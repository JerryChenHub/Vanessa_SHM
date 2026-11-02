import tempfile
from pathlib import Path
import unittest

from core.bdf_damage import read_cards, write_damage_regions, nastran_float


def small(name, *values):
    return f'{name:<8}'+''.join(f'{str(v):>8}' for v in values)+'\n'


def large(name, *values):
    return f'{name:<8}'+''.join(f'{str(v):>16}' for v in values)+'\n'


class DamageBdfTest(unittest.TestCase):
    def test_all_shell_stiffnesses_and_large_field_beam_preserve_mass_and_geometry(self):
        source = (small('MAT1',1,'1000.','400.','.25','2.7')
                  +small('MAT1',2,'2000.','800.','.25','3.0')
                  +small('PSHELL',10,1,'.1',2,'1.',1,'.833333')
                  +small('', '', '', 2)
                  +large('PBAR*',20,2,'.5275','.12455')
                  +large('*','.749325','.001758')
                  +small('CQUAD4',100,10,1,2,3,4)
                  +large('CBAR*',101,20,1,2)
                  +large('*','0.','1.','0.')
                  +small('GRID',1,'','0.','0.','0.'))
        with tempfile.TemporaryDirectory() as directory:
            base, out = Path(directory)/'base.bdf', Path(directory)/'case.bdf'
            base.write_text(source)
            _, changed = write_damage_regions([{1,2,3,4}],out,.5,base)
            self.assertEqual(changed,[100,101])
            self.assertEqual(base.read_text(),source)
            cards = read_cards(out.read_text().splitlines(keepends=True))
            by_key = {(c.name,int(c.fields[1])):c for c in cards}
            shell = by_key[('PSHELL',int(by_key[('CQUAD4',100)].fields[2]))]
            for field, original in [(2,1),(4,2),(6,1),(11,2)]:
                self.assertNotEqual(int(shell.fields[field]),original)
                new_mat = by_key[('MAT1',int(shell.fields[field]))]
                old_mat = by_key[('MAT1',original)]
                self.assertEqual(nastran_float(new_mat.fields[2]),.5*nastran_float(old_mat.fields[2]))
                self.assertEqual(nastran_float(new_mat.fields[3]),.5*nastran_float(old_mat.fields[3]))
                self.assertEqual(new_mat.fields[4:],old_mat.fields[4:])
            beam = by_key[('CBAR',101)]
            self.assertEqual(beam.fields[3:],["1","2","0.","1.","0.",""])
            prop = by_key[('PBAR',int(beam.fields[2]))]
            self.assertEqual(prop.fields[3:],by_key[('PBAR',20)].fields[3:])

    def test_different_regions_are_not_unioned_and_missing_properties_fail(self):
        source = small('MAT1',1,'1000.','400.','.25','2.7')+small('PBAR',10,1,'1.')+small('CBAR',100,10,1,2)
        with tempfile.TemporaryDirectory() as directory:
            base, out = Path(directory)/'base.bdf',Path(directory)/'case.bdf'
            base.write_text(source)
            self.assertEqual(write_damage_regions([{1},{2}],out,.5,base)[1],[])
            base.write_text(source.replace(small('PBAR',10,1,'1.'),''))
            with self.assertRaises(ValueError):
                write_damage_regions([{1,2}],out,.5,base)


if __name__ == '__main__':
    unittest.main()
