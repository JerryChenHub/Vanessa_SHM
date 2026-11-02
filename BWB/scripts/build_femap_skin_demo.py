"""Build both BWB skin demonstration models with SolidWorks-exported geometry.

No solver is run. Each part is a separate Femap database, in m-N-kg-s units.
The outer CAD skin is used as a zero-offset shell reference approximation.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from femap_connection import ROOT, connect, check, ids
from structured_skin_grid import generate, put_in_femap

OUT = ROOT / 'BWB/models/skin_structured'
CAD = ROOT / 'BWB/models/skin_demo/cad'
LOG = ROOT / 'workspace/BWB/femap_demo'

def all_ids(app, kind):
    s = app.feSet
    check(s.AddAll(kind), f'select type {kind}')
    return s, ids(s)

def bounds(app):
    solids, solid_ids = all_ids(app, 39)
    boxes = []
    body = app.feSolid
    for sid in solid_ids:
        check(body.Get(sid), 'solid')
        rc, box = body.BoundingBox()
        check(rc, 'solid bbox')
        boxes.append(box)
    a = np.asarray(boxes)
    return solids, np.r_[a[:, :3].min(0), a[:, 3:].max(0)]

def build(app, tag, size, chord_divisions=24):
    folder = OUT / f'part_{tag}'
    folder.mkdir(parents=True, exist_ok=True)
    check(app.feFileNew(), 'new model')
    geometry_import = app.feGeometryInterface
    geometry_import.ImportGeometryScale = 1.0  # Parasolid is in metres.
    geometry_import.ImportShowDialog = False
    check(geometry_import.Import(30, str(CAD/f'{tag}.x_t')), 'import geometry in metres')
    solids, box = bounds(app)
    # Independently measured exact SolidWorks z-extent, in metres.
    expected_span = 0.8730989581765075
    extent = float(box[5]-box[2])
    candidates = (1.0, 0.0254, 0.001)
    scale = min(candidates, key=lambda f: abs(extent*f-expected_span))
    if abs(extent*scale-expected_span) > 1e-4:
        raise RuntimeError(f'Unrecognized imported geometry scale: {box}')
    before = box.tolist()
    if scale != 1.0:
        check(app.feScale(39, solids.ID, (0.0, 0.0, 0.0), 0, scale, scale, scale), 'scale to metres')
    _, box = bounds(app)
    root_z = 0.0 if tag == '1' else 1.56
    all_surfaces, sids = all_ids(app, 5)
    skin = app.feSet
    surface = app.feSurface
    caps, surface_info = [], []
    for sid in sids:
        check(surface.Get(sid), 'get surface')
        rc, sb = surface.BoundingBox()
        check(rc, 'surface bbox')
        plane = surface.Planar() if surface.IsPlane() else None
        cap = bool(plane and plane[0] == -1 and abs(plane[1][2]-root_z) < 1e-6 and abs(plane[2][2]) > 0.999)
        surface_info.append(dict(id=sid,bounds=sb,plane=plane,excluded_root_cap=cap))
        if cap:
            caps.append(sid)
        else:
            check(skin.Add(sid), 'select skin surface')
    if len(caps) != 1:
        raise RuntimeError(f'Expected exactly one root cap: {caps}')
    print(f'Part {tag}: scale={scale}, bbox={box}, excluded root cap={caps}', flush=True)

    mat = app.feMatl
    mat.type = 0
    mat.title = 'DEMO aluminum: E69GPa nu0.33 rho2700kg/m3'
    mat.Ex, mat.Gx, mat.Nuxy, mat.Density = 69e9, 69e9/(2*1.33), .33, 2700.0
    check(mat.Put(1), 'material')
    prop = app.feProp
    prop.type, prop.matlID = 17, 1
    prop.title = 'DEMO skin: t=1mm; outer CAD reference; zero offset'
    prop.Setpval(0, .001)
    prop.SetExtraMatlID(0, 1)
    prop.SetExtraMatlID(1, 1)
    prop.Setpval(10, 1.0)
    prop.Setpval(11, 5/6)
    check(prop.Put(1), 'shell property')
    print(f'Part {tag}: building explicit span/chord grids, {chord_divisions} chord columns', flush=True)
    mesh_nodes, mesh_elements, layout = generate(app,size,chord_divisions)
    surface_element_map=put_in_femap(app,mesh_nodes,mesh_elements,layout['node_geometry'])
    layout['element_patches']={str(e['id']):e['patch'] for e in mesh_elements}
    layout['element_surface_ids']={str(e['id']):e['surface_id'] for e in mesh_elements}
    (folder/'mesh_layout.json').write_text(json.dumps(layout,indent=2),encoding='utf-8')
    rc, nn, nids, xyz = app.feNode.GetCoordArray(0)
    check(rc, 'mesh nodes')
    coords = np.asarray(xyz).reshape(-1,3)
    root_ids = [int(nids[i]) for i in np.flatnonzero(abs(coords[:,2]-root_z)<1e-7)]
    if len(root_ids)<10:
        raise RuntimeError(f'Insufficient root nodes: {len(root_ids)}')
    rootset = app.feSet
    check(rootset.AddArray(len(root_ids), tuple(root_ids)), 'root node selection')
    bcset = app.feBCSet
    bcset.title = f'DEMO fixed root z={root_z}m; DOF123456'
    check(bcset.Put(1), 'constraint set')
    bcset.Active = 1
    bc = app.feBCNode
    bc.SetID = 1
    check(bc.Add(rootset.ID, True, True, True, True, True, True), 'root clamp')
    am = app.feAnalysisMgr
    check(am.InitAnalysisMgr(36, 1, True), 'initialize NX static export settings')
    am.title = 'BWB shell mesh demonstration - NOT SOLVED'
    am.SetBCSet(0, 1)
    am.NasBulkWriteAllLoadBCSets = True
    check(am.Put(1), 'analysis export settings')
    am.Active = 1

    elements, eids = all_ids(app, 8)
    surface_counts = {}
    for sid in sids:
        check(surface.Get(sid), 'get surface mesh coverage')
        selected_elements = app.feSet
        surface.ElementsAsSet(selected_elements.ID, True)
        surface_counts[sid] = selected_elements.Count()
    quality = app.feElementQuality
    quality.GetModelDefaults()
    failed = quality.CheckNastranQuality(elements.ID, False, False)
    failures = ids(failed) if failed is not None else []
    aspect = []
    for eid in eids:
        rc, value = quality.GetAspectRatio(eid)
        if rc == -1:
            aspect.append((eid, value))

    rc, vid = app.feAppGetActiveView()
    check(rc, 'active view')
    view = app.feView
    check(view.Get(vid), 'get view')
    view.Mode, view.Contour, view.Deformed, view.BC = 0,0,0,1
    for idx in (9,20,25,26):
        view.SetDraw(idx,True)
    for idx in (2,3,4,5):
        view.SetDraw(idx,False)
    view.SetLabel(9,0)
    view.SetLabel(20,0)
    view.SetDraw(8,False)
    check(view.Put(vid),'mesh view')
    orientation=app.feViewOrient
    check(orientation.Get(vid),'get view orientation')
    orientation.Perspective=False
    orientation.SetRotationAngles(-65,20,-25)
    check(orientation.Put(vid),'view orientation')
    app.feAppVisible(True)
    app.feViewAutoscaleAll(vid,True)
    app.feViewRegenerate(vid)

    check(app.feFileSaveAs(False,str(folder/f'part_{tag}_skin.modfem')),'save Femap model')
    check(app.feFileWriteNastran(8,str(folder/f'part_{tag}_skin.bdf')),'export structural Nastran deck')
    report=dict(part=tag,units='m-N-kg-s',imported_bbox=before,import_scale_to_m=scale,
                geometry_bbox_m=box.tolist(),root_plane_z_m=root_z,excluded_cap_surfaces=caps,
                selected_skin_surface_ids=ids(skin),surface_inventory=surface_info,
                surface_element_counts={sid:len(surface_element_map.get(sid,[])) for sid in sids},
                femap_elements_all_nodes_on_surface=surface_counts,
                unmeshed_selected_surface_ids=[i for i in ids(skin) if not surface_element_map.get(i)],
                target_size_m=size,nodes=nn,elements=len(eids),root_nodes=root_ids,
                mesh_method=layout['method'],chord_divisions=chord_divisions,
                span_divisions=layout['span_divisions'],
                node_merge_tolerance_m=1e-8,node_merge_tolerance_auto=False,
                femap_nastran_quality_failed_elements=failures,
                femap_aspect_ratio_max=max(aspect,key=lambda x:x[1]) if aspect else None,
                femap_aspect_ratio_p95=float(np.percentile([a[1] for a in aspect],95)) if aspect else None,
                material=dict(E_Pa=69e9,nu=.33,rho_kg_m3=2700),thickness_m=.001,
                shell_reference='Outer CAD surface with zero offset; demonstration approximation',
                solved=False,includes_aerodynamics=False)
    (folder/'build_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('surface_inventory','root_nodes')}),flush=True)
    return report

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--parts',nargs='+',default=['1','2'],choices=['1','2'])
    parser.add_argument('--size',type=float,default=.020)
    parser.add_argument('--chord-divisions',type=int,default=24)
    args=parser.parse_args()
    if args.size<=0 or args.chord_divisions<4:
        parser.error('--size must be positive and --chord-divisions must be at least 4')
    app=connect(True)
    previous=(app.DialogAutoAnswer,app.DialogAutoSkip,app.DialogAutoSkipMsg)
    try:
        app.DialogAutoAnswer=app.DialogAutoSkip=app.DialogAutoSkipMsg=1
        app.Tool_MergeToleranceAuto=False
        app.Tool_MergeTolerance=1e-8
        for tag in args.parts:
            build(app,tag,args.size,args.chord_divisions)
    finally:
        app.DialogAutoAnswer,app.DialogAutoSkip,app.DialogAutoSkipMsg=previous

if __name__=='__main__':
    main()
