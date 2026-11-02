"""Fixed/large-field BDF damage editing without rewriting unchanged geometry."""
from dataclasses import dataclass
import math
from pathlib import Path
import re
from .paths import require_working_output

ELEMENT_NODES = {'CBAR':2, 'CBEAM':2, 'CROD':2, 'CTRIA3':3, 'CQUAD4':4, 'CTETRA':4, 'CHEXA':8}


def nastran_float(text):
    text = str(text).strip().replace('D','E').replace('d','e')
    if not text:
        return None
    if 'e' not in text.lower():
        text = re.sub(r'(?<=[\d.])([+-]\d+)$', r'e\1', text)
    return float(text)


def real8(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('Nonfinite BDF real')
    if value == 0:
        return '0.'
    candidates = []
    for precision in range(8):
        text = f'{value:.{precision}f}'
        if '.' in text:
            text = text.rstrip('0').rstrip('.')
        if '.' not in text:
            text += '.'
        candidates.append(text)
        mantissa, exponent = f'{value:.{precision}E}'.split('E')
        mantissa = mantissa.rstrip('0').rstrip('.') if '.' in mantissa else mantissa
        if '.' not in mantissa:
            mantissa += '.'
        candidates.append(f'{mantissa}{int(exponent):+d}')
    candidates = [s for s in candidates if len(s)<=8 and nastran_float(s)!=0]
    if not candidates:
        raise ValueError(f'Cannot represent {value} in an 8-column real field')
    return min(candidates, key=lambda s:abs(nastran_float(s)-value))


@dataclass
class Card:
    name: str
    fields: list
    locations: list
    line_numbers: list

    def edited_lines(self, lines, changes):
        result = {i:lines[i] for i in self.line_numbers}
        for field, value in changes.items():
            line_index, offset, width = self.locations[field]
            text = str(value)
            if len(text)>width:
                raise ValueError(f'{self.name} field {field} does not fit {width} columns')
            original = result[line_index].rstrip('\r\n').ljust(offset+width)
            result[line_index] = original[:offset]+text.rjust(width)+original[offset+width:]+'\n'
        return result


def read_cards(lines):
    cards, current, width = [], None, 8
    for index, line in enumerate(lines):
        if not line.strip() or line.lstrip().startswith('$'):
            continue
        name = line[:8].strip()
        continuation = not name or name.startswith(('+','*'))
        if not continuation:
            if ',' in line[:8]:
                raise ValueError('This mesh reader requires fixed or large field BDF records')
            width = 16 if name.endswith('*') else 8
            current = Card(name.rstrip('*'), [name.rstrip('*')], [None], [])
            cards.append(current)
        elif current is None:
            raise ValueError('Continuation without a card')
        current.line_numbers.append(index)
        for offset in range(8,72,width):
            current.fields.append(line[offset:offset+width].strip())
            current.locations.append((index,offset,width))
    return cards


def write_damage_regions(selected_nodes_groups, output_path, damage_ratio, base_bdf_path):
    source, target = Path(base_bdf_path), require_working_output(output_path)
    if source.resolve() == target.resolve():
        raise ValueError('The damaged case must not overwrite its base mesh')
    if not 0 < damage_ratio <= 1:
        raise ValueError('damage_ratio is a remaining stiffness scale in (0,1]')
    lines = source.read_text().splitlines(keepends=True)
    cards = read_cards(lines)
    groups = [set(map(int,g)) for g in selected_nodes_groups]
    materials = {int(c.fields[1]):c for c in cards if c.name=='MAT1'}
    properties = {int(c.fields[1]):c for c in cards if c.name in {'PSHELL','PBAR'}}
    selected = []
    for card in cards:
        if card.name in ELEMENT_NODES:
            nodes = [int(n) for n in card.fields[3:3+ELEMENT_NODES[card.name]] if n]
            if nodes and any(all(n in group for n in nodes) for group in groups):
                selected.append(card)
    for card in selected:
        pid = int(card.fields[2])
        if pid not in properties:
            raise ValueError(f'Element {card.fields[1]} uses unsupported/missing property {pid}')
        prop = properties[pid]
        indices = (2,4,6,11) if prop.name=='PSHELL' else (2,)
        for field in indices:
            if field<len(prop.fields) and prop.fields[field] and int(prop.fields[field])>0 and int(prop.fields[field]) not in materials:
                raise ValueError(f'Property {pid} uses unsupported/missing MAT1')
    next_pid, next_mid = max(properties)+1000, max(materials)+1000
    mapping, additions = {}, []
    for pid in sorted({int(c.fields[2]) for c in selected}):
        prop = properties[pid]
        additions.append('$ Damaged material/property\n')
        indices = (2,4,6,11) if prop.name=='PSHELL' else (2,)
        property_edits, material_map = {1:next_pid}, {}
        for field in indices:
            if field>=len(prop.fields) or not prop.fields[field] or int(prop.fields[field])<=0:
                continue
            old_mid = int(prop.fields[field])
            if old_mid not in material_map:
                material = materials[old_mid]
                edits = {1:next_mid}
                for modulus_field in (2,3):
                    value = nastran_float(material.fields[modulus_field])
                    if value is not None:
                        width = material.locations[modulus_field][2]
                        edits[modulus_field] = real8(value*damage_ratio) if width==8 else f'{value*damage_ratio:.9E}'
                additions.extend(material.edited_lines(lines,edits).values())
                material_map[old_mid] = next_mid
                next_mid += 1
            property_edits[field] = material_map[old_mid]
        additions.extend(prop.edited_lines(lines,property_edits).values())
        mapping[pid] = next_pid
        next_pid += 1
    output = lines.copy()
    for card in selected:
        for index, line in card.edited_lines(lines,{2:mapping[int(card.fields[2])]}).items():
            output[index] = line
    insertion = next((c.line_numbers[0] for c in cards if c.name=='GRID'),len(output))
    output[insertion:insertion] = additions
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(''.join(output))
    print('modified elements:',len(selected))
    print('new properties:',len(mapping))
    print('pid map:',mapping)
    print('saved:',target)
    return str(target), [int(c.fields[1]) for c in selected]
