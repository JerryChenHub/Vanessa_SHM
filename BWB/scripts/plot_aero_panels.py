"""Render actual CAERO1 cards beside their source skin outlines, without solving."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

try:
    from .check_aero_panels import grouped_cards, read_aero_geometry
    from .inspect_shell_bdf import load_model
except ImportError:
    from check_aero_panels import grouped_cards, read_aero_geometry
    from inspect_shell_bdf import load_model

ROOT = Path(__file__).resolve().parents[2]
DEFAULT = ROOT/'BWB/models/aero_panels'
COLORS = {'main':'#559bcf','transition':'#e9b253','bend':'#a382cb','winglet':'#55ad97'}
LABELS = {'main':'Main wing','transition':'Transition','bend':'Curved bend','winglet':'Winglet'}


def plotted(points):
    return np.asarray(points)[..., [2,0,1]]*1000


def render(folders, output, support_points=False):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':'#243c4f',
                         'axes.labelcolor':'#506b7e','xtick.color':'#617b8d','ytick.color':'#617b8d'})
    fig=plt.figure(figsize=(16,8.8),facecolor='#f7f9fc')
    fig.text(.045,.947,'BWB | '+('Spline support locations' if support_points else 'Generated aerodynamic panels'),
             fontsize=23,weight='bold')
    fig.text(.045,.901,'Two independent components  |  Native CAERO1 geometry  |  Source skin boundaries shown in grey',
             fontsize=12,color='#536d7e')
    for index, folder in enumerate(map(Path, folders)):
        manifest=json.loads((folder/'aero_manifest.json').read_text(encoding='utf-8'))
        geometry=read_aero_geometry(grouped_cards(folder/manifest['files']['aero']))
        bdf_source=ROOT/next(p['path'] for p in manifest['source_files'] if p['path'].endswith('.bdf'))
        layout_source=ROOT/next(p['path'] for p in manifest['source_files'] if p['path'].endswith('mesh_layout.json'))
        skin=load_model(bdf_source)
        layout=json.loads(layout_source.read_text(encoding='utf-8'))
        coordinates={nid:np.asarray(n['xyz']) for nid,n in skin['nodes'].items()}
        ax=fig.add_subplot(1,2,index+1,projection='3d',facecolor='#f7f9fc')
        for zone in COLORS:
            for side in ['lower','upper']:
                grid=np.array([[coordinates[n] for n in row] for row in layout['patch_node_ids'][zone+'_'+side]])
                for edge in [grid[:,0],grid[:,-1],grid[0],grid[-1]]:
                    p=plotted(edge)
                    ax.plot(*p.T,color='#9ba9b4',linewidth=.7,alpha=.7)
        for panel in manifest['panels']:
            macro=geometry['macros'][panel['eid']]
            polys=[plotted(geometry['boxes'][i]['corners']) for i in macro['ids']]
            ax.add_collection3d(Poly3DCollection(polys,facecolors=COLORS[panel['zone']],
                edgecolors='#304d61',linewidths=.35,alpha=.26 if support_points else .85))
            g=macro['points']
            for edge in [g[0],g[-1],g[:,0],g[:,-1]]:
                p=plotted(edge)
                ax.plot(*p.T,color='#263e52',linewidth=.8,alpha=.7)
        if support_points:
            selected=sorted({n for p in manifest['panels'] for n in p['structural_node_ids']})
            p=plotted([coordinates[n] for n in selected])
            ax.scatter(*p.T,c='#c73356',s=10,depthshade=False,edgecolors='white',linewidths=.3)
        root=plotted([coordinates[n] for n in layout['patch_node_ids']['main_upper'][0]])
        ax.plot(*root.T,color='#cd4350',linewidth=2.8)
        p=plotted(list(coordinates.values()))
        low,high=p.min(axis=0),p.max(axis=0)
        size=high-low
        for k,setter in enumerate([ax.set_xlim,ax.set_ylim,ax.set_zlim]):
            margin=max(size[k]*.06,size.max()*.02)
            setter(low[k]-margin,high[k]+margin)
        ax.set_box_aspect(np.maximum(size,size.max()*.06))
        ax.set_proj_type('ortho')
        ax.view_init(elev=27,azim=-65 if index==0 else -115)
        ax.set_xlabel('CAD Z / span (mm)',labelpad=10)
        ax.set_ylabel('CAD X / chord (mm)',labelpad=10)
        ax.set_zlabel('CAD Y / height (mm)',labelpad=8)
        ax.tick_params(labelsize=8,pad=1)
        ax.set_title(f"PART {manifest['part']}  |  {len(geometry['macros'])} regions  |  {len(geometry['boxes'])} boxes",
                     loc='left',fontsize=12,y=.94)
        for pane in [ax.xaxis.pane,ax.yaxis.pane,ax.zaxis.pane]:
            pane.set_facecolor((1,1,1,.2))
            pane.set_edgecolor((.8,.85,.9,.4))
        for a in [ax.xaxis,ax.yaxis,ax.zaxis]:
            a._axinfo['grid'].update(color=(.75,.8,.85,.35),linewidth=.6)
        x=.05 if index==0 else .54
        fig.text(x,.175,'10 chordwise boxes; span divisions: 39 + 3 + 8 + 9',fontsize=10,color='#536b7c')
        fig.text(x,.145,f"Source: part_{manifest['part']}/aero_cards.bdf",fontsize=9,color='#718495')
    handles=[Patch(facecolor=COLORS[z],edgecolor='#304d61',label=LABELS[z]) for z in COLORS]
    handles.append(Line2D([0],[0],color='#cd4350',linewidth=3,label='Fixed structural root'))
    if support_points:
        handles.append(Line2D([0],[0],marker='o',linestyle='',color='#c73356',label='Spline GRID points'))
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.04,.075),ncol=len(handles),frameon=False,fontsize=9)
    fig.text(.045,.055,'Actual exported CAERO1 cards, independently expanded for display. Original CAD coordinates; these are not an assembly.',
             fontsize=9,color='#718495')
    fig.text(.045,.025,'W2GJ mean-camber corrections are stored separately in dmicards.bdf. No pressure, lift, strain or SOL144 result is shown.',
             fontsize=9,color='#718495')
    fig.subplots_adjust(left=.015,right=.985,top=.86,bottom=.21,wspace=.04)
    output=Path(output)
    output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(output,dpi=175,facecolor=fig.get_facecolor())
    plt.close(fig)
    print(output)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',type=Path,default=DEFAULT)
    args=p.parse_args()
    folders=[args.root/f'part_{part}' for part in (1,2)]
    render(folders,args.root/'aero_panel_overview.png')
    render(folders,args.root/'aero_spline_support.png',support_points=True)


if __name__=='__main__':
    main()
