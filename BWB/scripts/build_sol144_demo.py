"""Build reproducible, independently clamped BWB SOL144 demonstration cases."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil

from .build_bwb_aero_panels import ROOT, provenance, real, write_card

DEFAULT_OUTPUT = ROOT / 'BWB/cases/sol144_demo'
INPUT_NAMES = ('structural_cards.bdf', 'aero_cards.bdf', 'spline_cards.bdf', 'dmicards.bdf')


def build(output=DEFAULT_OUTPUT, velocity=20., density=1.225, sound_speed=340.3, alpha_deg=2.):
    output = Path(output).resolve()
    if not output.is_relative_to(ROOT / 'BWB') and not output.is_relative_to(ROOT / 'workspace/BWB'):
        raise ValueError('Output must be owned by BWB')
    if not all(math.isfinite(x) for x in (velocity, density, sound_speed, alpha_deg)):
        raise ValueError('Flight parameters must be finite')
    if min(velocity, density, sound_speed) <= 0 or velocity >= sound_speed:
        raise ValueError('This demonstration requires positive, subsonic flight parameters')
    q = .5 * density * velocity**2
    mach = velocity / sound_speed
    for part in (1, 2):
        source = ROOT / f'BWB/models/aero_panels/part_{part}'
        dest = output / f'part_{part}'
        dest.mkdir(parents=True, exist_ok=True)
        for name in INPUT_NAMES:
            shutil.copyfile(source / name, dest / name)
        condition = '$ Assumed aerodynamic-only fixed-root demonstration; SI units.\n'
        condition += write_card('AESTAT', [1, 'ANGLEA'])
        condition += write_card('AESTAT', [2, 'SIDES'])
        # On this Simcenter release a blank second label/value pair caused
        # AEQR to be ignored. Explicit zero sideslip preserves the flight case
        # and demonstrably activates the aeroelastic stiffness contribution.
        condition += write_card('TRIM', [1, real(mach), real(q), 'ANGLEA', real(math.radians(alpha_deg)), 'SIDES', real(0.), real(1.)])
        (dest / 'flight_condition.bdf').write_text(condition, encoding='ascii')
        deck = f'''$ BWB part {part}: independent clamped component, not full-aircraft trim.
SOL 144
CEND
TITLE = BWB PART {part} ASSUMED LOW SPEED AEROELASTIC STRAIN DEMO
ECHO = NONE
SUBCASE 1
  LABEL = ASSUMED_LOW_SPEED_ELASTIC
  SPC = 1
  TRIM = 1
  DISPLACEMENT(PRINT,PLOT,SORT1,REAL) = ALL
  SPCFORCES(PRINT,PLOT,SORT1,REAL) = ALL
  STRESS(PRINT,PLOT,SORT1,REAL,VONMISES,CENTER) = ALL
  STRAIN(PRINT,PLOT,SORT1,REAL,FIBER,CENTER) = ALL
  APRES = ALL
  AEROF = ALL
BEGIN BULK
INCLUDE 'structural_cards.bdf'
INCLUDE 'aero_cards.bdf'
INCLUDE 'spline_cards.bdf'
INCLUDE 'dmicards.bdf'
INCLUDE 'flight_condition.bdf'
ENDDATA
'''
        (dest / 'bwb_sol144.bdf').write_text(deck, encoding='ascii')
        manifest = {
            'schema_version': 1, 'part': part, 'units': 'm-N-kg-s',
            'flight': {'velocity_m_s': velocity, 'density_kg_m3': density,
                       'sound_speed_m_s': sound_speed, 'mach': mach,
                       'dynamic_pressure_Pa': q, 'angle_of_attack_deg': alpha_deg,
                       'angle_of_attack_rad': math.radians(alpha_deg), 'sideslip_deg': 0.},
            'assumptions': [
                'Each part is an independent root-clamped component, with no symmetry image.',
                'Linear static aeroelasticity (TRIM AEQR=1); camber/twist W2GJ retained.',
                'Sea-level-like density and sound speed are assumed constants, not site measurements.',
                'Aluminium E=69 GPa, nu=0.33, rho=2700 kg/m^3; uniform 1 mm skin.',
                '49 root nodes constrained in all six DOF; no ribs, spars, attached masses or controls.',
                'Aerodynamic loading only; gravity, maneuver acceleration and residual stress omitted.',
                'No prescribed CL: Nastran solves panel pressure and coupled elastic displacement.',
                'No free-flight force/moment equilibrium; root reactions support the aerodynamic load.',
                'Strains are linear shell centroid fiber strains in element coordinates.',
            ],
            'source_inputs': [provenance(source / name) for name in INPUT_NAMES],
            'generator': provenance(Path(__file__)),
            'nasa_read_only_reference': [provenance(ROOT / f'NASA_SHM/core/base/{name}')
                                       for name in ('sol144.bdf', 'subcases.inc', 'trim_cards.bdf')],
            'nasa_adaptation': 'SOL144, SPC/TRIM and shell output organization only. NASA inch/weight unit parameters and numeric trim values are not copied.',
            'inputs': [provenance(dest / name) for name in (*INPUT_NAMES, 'flight_condition.bdf', 'bwb_sol144.bdf')],
        }
        (dest / 'case_manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--velocity', type=float, default=20.)
    parser.add_argument('--density', type=float, default=1.225)
    parser.add_argument('--sound-speed', type=float, default=340.3)
    parser.add_argument('--alpha-deg', type=float, default=2.)
    args = parser.parse_args()
    print(build(args.output, args.velocity, args.density, args.sound_speed, args.alpha_deg))
