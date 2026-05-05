import os
import csv
bdf_path = r"/base_file/FEM_only.bdf"

output_dir = os.path.dirname(bdf_path)
nodes_csv = os.path.join(output_dir, "nodes.csv")
elements_csv = os.path.join(output_dir, "elements.csv")


def to_float(s):
    s = s.strip()
    if not s:
        return None

    if "e" not in s.lower():
        for i in range(1, len(s)):
            if s[i] in "+-" and s[i - 1].isdigit():
                s = s[:i] + "e" + s[i:]
                break

    return float(s)


def fixed_fields(line, width=8):
    return [line[i:i + width].strip() for i in range(0, len(line.rstrip("\n")), width)]


def parse_bdf(file_path):
    nodes = []
    elements = []

    element_cards = {
        "CBAR": 2,
        "CBEAM": 2,
        "CROD": 2,
        "CTRIA3": 3,
        "CQUAD4": 4,
        "CTETRA": 4,
        "CHEXA": 8,
    }

    with open(file_path, "r", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("$"):
                continue

            card = line[:8].strip()

            if card in ("", "*", "+"):
                continue

            fields = fixed_fields(line)
            card = fields[0]

            if card == "GRID":
                try:
                    node_id = int(fields[1])
                    x = to_float(fields[3])
                    y = to_float(fields[4])
                    z = to_float(fields[5])
                    nodes.append([node_id, x, y, z])
                except Exception:
                    pass

            elif card in element_cards:
                try:
                    elem_id = int(fields[1])
                    prop_id = int(fields[2])
                    n_nodes = element_cards[card]

                    conn = []
                    for i in range(n_nodes):
                        val = fields[3 + i]
                        if val:
                            conn.append(int(val))

                    elements.append([elem_id, card, prop_id] + conn)
                except Exception:
                    pass

    return nodes, elements


def write_nodes(nodes, filepath):
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["node_id", "x", "y", "z"])
        writer.writerows(nodes)


def write_elements(elements, filepath):
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)

        if not elements:
            writer.writerow(["elem_id", "type", "prop_id"])
            return

        max_conn = max(len(e) - 3 for e in elements)
        header = ["elem_id", "type", "prop_id"] + [f"n{i+1}" for i in range(max_conn)]
        writer.writerow(header)

        for e in elements:
            row = e + [""] * (3 + max_conn - len(e))
            writer.writerow(row)

def fixed_fields(line, width=8):
    return [line[i:i + width].strip() for i in range(0, len(line.rstrip("\n")), width)]


def expand_thru(tokens):
    nodes = []
    i = 0

    while i < len(tokens):
        token = tokens[i].upper()

        if i + 2 < len(tokens) and tokens[i + 1].upper() == "THRU":
            start = int(tokens[i])
            end = int(tokens[i + 2])
            nodes.extend(range(start, end + 1))
            i += 3
        else:
            try:
                nodes.append(int(tokens[i]))
            except Exception:
                pass
            i += 1

    return nodes


def parse_bdf_cards(file_path):
    cards = []
    current = None

    with open(file_path, "r", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("$"):
                continue

            fields = fixed_fields(line)
            card = fields[0]

            if card:
                if current is not None:
                    cards.append(current)
                current = fields

            else:
                if current is not None:
                    current.extend(fields[1:])

        if current is not None:
            cards.append(current)

    return cards


def parse_fixed_nodes(file_path):
    fixed_nodes = set()

    cards = parse_bdf_cards(file_path)

    for fields in cards:
        card = fields[0]

        if card == "SPC":
            # SPC: SID, G1, C1, D1, G2, C2, D2...
            i = 2
            while i < len(fields):
                try:
                    if fields[i]:
                        fixed_nodes.add(int(fields[i]))
                except Exception:
                    pass
                i += 3

        elif card == "SPC1":
            # SPC1: SID, C, G1, G2, ...
            node_tokens = [x for x in fields[3:] if x]
            fixed_nodes.update(expand_thru(node_tokens))

    return sorted(fixed_nodes)


def write_fixed_nodes(fixed_nodes, filepath):
    with open(filepath, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["node_id"])
        for node_id in fixed_nodes:
            writer.writerow([node_id])


import os
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


def plot_bdf_mesh_3d(data_dir="data", show_elements=True, show_faces=True):
    nodes_path = os.path.join(data_dir, "nodes.csv")
    elements_path = os.path.join(data_dir, "elements.csv")
    fixed_path = os.path.join(data_dir, "fixed_nodes.csv")

    nodes_df = pd.read_csv(nodes_path)
    elements_df = pd.read_csv(elements_path)
    fixed_df = pd.read_csv(fixed_path)

    # node_id -> xyz
    node_xyz = {
        int(row["node_id"]): (row["x"], row["y"], row["z"])
        for _, row in nodes_df.iterrows()
    }

    fixed_nodes = set(fixed_df["node_id"].astype(int).tolist())

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection="3d")

    ax.scatter(
        nodes_df["x"],
        nodes_df["y"],
        nodes_df["z"],
        s=2,
        alpha=0.35,
        label="All nodes"
    )

    fixed_plot = nodes_df[nodes_df["node_id"].astype(int).isin(fixed_nodes)]

    ax.scatter(
        fixed_plot["x"],
        fixed_plot["y"],
        fixed_plot["z"],
        s=35,
        marker="^",
        c="black",
        label="Fixed nodes"
    )

    # ===== 组 element =====
    faces = []

    node_cols = [c for c in elements_df.columns if c.startswith("n")]

    for _, elem in elements_df.iterrows():
        elem_type = str(elem["type"]).strip()

        conn = []
        for col in node_cols:
            if pd.notna(elem[col]):
                try:
                    nid = int(elem[col])
                    if nid in node_xyz:
                        conn.append(nid)
                except Exception:
                    pass

        if len(conn) < 2:
            continue

        pts = [node_xyz[nid] for nid in conn]

        if elem_type in ["CBAR", "CBEAM", "CROD"]:
            if show_elements and len(pts) >= 2:
                xs, ys, zs = zip(*pts[:2])
                ax.plot(xs, ys, zs, linewidth=0.35, alpha=0.45)

        elif elem_type in ["CTRIA3", "CQUAD4"]:
            if show_elements:
                closed_pts = pts + [pts[0]]
                xs, ys, zs = zip(*closed_pts)
                ax.plot(xs, ys, zs, linewidth=0.25, alpha=0.35)

            if show_faces:
                faces.append(pts)

        else:
            if show_elements:
                closed_pts = pts + [pts[0]]
                xs, ys, zs = zip(*closed_pts)
                ax.plot(xs, ys, zs, linewidth=0.25, alpha=0.25)

    # ===== 画面片 =====
    if show_faces and faces:
        poly = Poly3DCollection(
            faces,
            alpha=0.15,
            edgecolor="none"
        )
        ax.add_collection3d(poly)

    set_axes_equal(ax)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("BDF Mesh with Fixed Nodes")
    ax.legend()

    plt.tight_layout()
    plt.show()


def set_axes_equal(ax):
    """
    """
    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    y_range = abs(y_limits[1] - y_limits[0])
    z_range = abs(z_limits[1] - z_limits[0])

    max_range = max(x_range, y_range, z_range) / 2.0

    x_middle = sum(x_limits) / 2.0
    y_middle = sum(y_limits) / 2.0
    z_middle = sum(z_limits) / 2.0

    ax.set_xlim3d([x_middle - max_range, x_middle + max_range])
    ax.set_ylim3d([y_middle - max_range, y_middle + max_range])
    ax.set_zlim3d([z_middle - max_range, z_middle + max_range])

def print_bounds(data_dir="data"):
    import os, pandas as pd
    n=pd.read_csv(os.path.join(data_dir,"nodes.csv"))
    f=pd.read_csv(os.path.join(data_dir,"fixed_nodes.csv"))
    fn=n[n["node_id"].isin(f["node_id"])]
    def b(df): return df[["x","y","z"]].min().values, df[["x","y","z"]].max().values
    nmin,nmax=b(n); fmin,fmax=b(fn)
    print("ALL min:",nmin,"max:",nmax)
    print("FIX min:",fmin,"max:",fmax)



def plot_nodes_fixed_box(data_dir="data"):
    import os, pandas as pd, matplotlib.pyplot as plt
    from itertools import product, combinations

    n=pd.read_csv(os.path.join(data_dir,"nodes.csv"))
    e=pd.read_csv(os.path.join(data_dir,"elements.csv"))
    f=pd.read_csv(os.path.join(data_dir,"fixed_nodes.csv"))

    fn=n[n["node_id"].isin(f["node_id"])]

    node_cols=[c for c in e.columns if c.startswith("n")]
    used=set()
    for col in node_cols:
        used.update(e[col].dropna().astype(int).tolist())

    floating=n[~n["node_id"].isin(used)]

    fig=plt.figure(figsize=(12,9))
    ax=fig.add_subplot(111,projection="3d")

    ax.scatter(n["x"],n["y"],n["z"],s=2,alpha=0.3)
    ax.scatter(fn["x"],fn["y"],fn["z"],s=35,c="black",marker="^")
    ax.scatter(floating["x"],floating["y"],floating["z"],s=40,facecolors='none',edgecolors='red',marker='o')

    mn=fn[["x","y","z"]].min()-10
    mx=fn[["x","y","z"]].max()+10

    corners=list(product([mn["x"],mx["x"]],[mn["y"],mx["y"]],[mn["z"],mx["z"]]))
    for a,b in combinations(corners,2):
        if sum(a[i]!=b[i] for i in range(3))==1:
            ax.plot([a[0],b[0]],[a[1],b[1]],[a[2],b[2]],linewidth=1.2)

    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    plt.tight_layout()
    plt.show()

def plot_pick_nodes_fixed_box(data_dir="data", offset=30):
    import os, pandas as pd, matplotlib.pyplot as plt
    from itertools import product, combinations

    n=pd.read_csv(os.path.join(data_dir,"nodes.csv"))
    e=pd.read_csv(os.path.join(data_dir,"elements.csv"))
    f=pd.read_csv(os.path.join(data_dir,"fixed_nodes.csv"))

    fn=n[n["node_id"].isin(f["node_id"])]

    node_cols=[c for c in e.columns if c.startswith("n")]
    used=set()
    for col in node_cols:
        used.update(e[col].dropna().astype(int).tolist())
    floating=n[~n["node_id"].isin(used)]

    fig=plt.figure(figsize=(12,9))
    ax=fig.add_subplot(111,projection="3d")

    sc=ax.scatter(n["x"],n["y"],n["z"],s=4,alpha=0.35,picker=True,pickradius=5)
    ax.scatter(fn["x"],fn["y"],fn["z"],s=45,c="black",marker="^")
    ax.scatter(floating["x"],floating["y"],floating["z"],s=50,facecolors="none",edgecolors="red",marker="o")

    mn=fn[["x","y","z"]].min()-offset
    mx=fn[["x","y","z"]].max()+offset
    corners=list(product([mn["x"],mx["x"]],[mn["y"],mx["y"]],[mn["z"],mx["z"]]))

    for a,b in combinations(corners,2):
        if sum(a[i]!=b[i] for i in range(3))==1:
            ax.plot([a[0],b[0]],[a[1],b[1]],[a[2],b[2]],linewidth=1.2)

    txt=ax.text2D(0.02,0.96,"Pick a node",transform=ax.transAxes)

    def on_pick(event):
        if event.artist!=sc or len(event.ind)==0:
            return
        i=event.ind[0]
        r=n.iloc[i]
        txt.set_text(f"id={int(r.node_id)}  x={r.x:.6g}  y={r.y:.6g}  z={r.z:.6g}")
        print(f"node_id={int(r.node_id)}, x={r.x}, y={r.y}, z={r.z}")
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect("pick_event",on_pick)

    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    plt.tight_layout()
    plt.show()


###从这里开始，开始选椭圆

def get_floating_nodes(data_dir="data"):
    import os, pandas as pd
    n=pd.read_csv(os.path.join(data_dir,"nodes.csv"))
    e=pd.read_csv(os.path.join(data_dir,"elements.csv"))
    node_cols=[c for c in e.columns if c.startswith("n")]
    used=set()
    for col in node_cols:
        used.update(e[col].dropna().astype(int).tolist())
    return set(n.loc[~n["node_id"].isin(used),"node_id"].astype(int))


def get_box_nodes(data_dir="data", offset=30):
    import os, pandas as pd
    n=pd.read_csv(os.path.join(data_dir,"nodes.csv"))
    f=pd.read_csv(os.path.join(data_dir,"fixed_nodes.csv"))
    fn=n[n["node_id"].isin(f["node_id"])]
    mn=fn[["x","y","z"]].min()-offset
    mx=fn[["x","y","z"]].max()+offset
    inside=(n["x"].between(mn["x"],mx["x"]) &
            n["y"].between(mn["y"],mx["y"]) &
            n["z"].between(mn["z"],mx["z"]))
    return set(n.loc[inside,"node_id"].astype(int))


def select_nodes_in_ellipsoid(node=1000,a=30,b=30,c=30,Ome=None,decay_rate=0.5,data_dir="data",box_offset=30):
    import os, numpy as np, pandas as pd

    if Ome is None:
        Ome=np.eye(3)

    n=pd.read_csv(os.path.join(data_dir,"nodes.csv"))
    center=n[n["node_id"]==node]
    if center.empty:
        raise ValueError(f"node {node} not found")

    p0=center[["x","y","z"]].iloc[0].to_numpy(float)
    P=n[["x","y","z"]].to_numpy(float)

    floating=get_floating_nodes(data_dir)
    box_nodes=get_box_nodes(data_dir,offset=box_offset)
    exclude=floating | box_nodes

    d=P-p0
    q=d @ np.asarray(Ome,float)
    val=(q[:,0]/a)**2+(q[:,1]/b)**2+(q[:,2]/c)**2
    weight=np.exp(-decay_rate*val)

    mask=(val<=1.0) & (~n["node_id"].astype(int).isin(exclude))
    selected=n.loc[mask,"node_id"].astype(int).tolist()

    return selected

def plot_selected_nodes(node=1000,a=30,b=30,c=30,Ome=None,decay_rate=0.5,data_dir="data",box_offset=30):
    import os, numpy as np, pandas as pd, matplotlib.pyplot as plt

    if Ome is None:
        Ome=np.eye(3)

    n=pd.read_csv(os.path.join(data_dir,"nodes.csv"))
    selected=select_nodes_in_ellipsoid(node,a,b,c,Ome,decay_rate,data_dir,box_offset)

    sel=n[n["node_id"].isin(selected)]
    center=n[n["node_id"]==node]

    fig=plt.figure(figsize=(10,8))
    ax=fig.add_subplot(111,projection="3d")

    ax.scatter(n["x"],n["y"],n["z"],s=2,alpha=0.2)
    ax.scatter(sel["x"],sel["y"],sel["z"],s=20,c="blue")
    ax.scatter(center["x"],center["y"],center["z"],s=80,c="green",marker="*")

    ax.set_xlabel("X"); ax.set_ylabel("Y"); ax.set_zlabel("Z")
    plt.tight_layout()
    plt.show()

def write_damaged_bdf(selected_nodes,bdf_in="data/FEM_only.bdf",damage_ratio=0.5,out_dir="data"):
    import os,time,re

    with open(bdf_in,"r",errors="ignore") as f:
        lines=f.readlines()

    def ff(line): return [line[i:i+8].strip() for i in range(0,len(line.rstrip("\n")),8)]
    def field(v): return str(v)[:8].rjust(8)
    def card(v): return str(v)[:8].ljust(8)
    def make(fields): return card(fields[0])+"".join(field(x) for x in fields[1:]).rstrip()+"\n"

    def to_float(s):
        s=str(s).strip().replace("D","E").replace("d","e")
        if not s: return None
        if "e" in s.lower(): return float(s)
        m=re.fullmatch(r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))([+-]\d+)",s)
        if m: return float(m.group(1)+"e"+m.group(2))
        return float(s)

    selected_nodes=set(map(int,selected_nodes))
    elem_cards={"CBAR":2,"CBEAM":2,"CROD":2,"CTRIA3":3,"CQUAD4":4,"CTETRA":4,"CHEXA":8}
    prop_cards={"PSHELL","PBAR"}

    props={}
    mats={}
    selected_elems=[]

    for line in lines:
        flds=ff(line)
        c=flds[0]

        if c in elem_cards:
            eid=int(flds[1])
            pid=int(flds[2])
            conn=[int(x) for x in flds[3:3+elem_cards[c]] if x]
            if conn and all(x in selected_nodes for x in conn):
                selected_elems.append((eid,pid))

        elif c in prop_cards:
            props[int(flds[1])]={"fields":flds.copy(),"mid":int(flds[2])}

        elif c=="MAT1":
            mats[int(flds[1])]=flds.copy()

    next_pid=max(props.keys())+1000
    next_mid=max(mats.keys())+1000

    pid_map={}
    add_lines=[]

    for old_pid in sorted(set(pid for _,pid in selected_elems)):
        if old_pid not in props: continue
        old_mid=props[old_pid]["mid"]
        if old_mid not in mats: continue

        new_pid=next_pid; next_pid+=1
        new_mid=next_mid; next_mid+=1

        mf=mats[old_mid].copy()
        while len(mf)<7: mf.append("")
        mf[1]=str(new_mid)
        mf[2]=f"{to_float(mf[2])*damage_ratio:.6g}"

        pf=props[old_pid]["fields"].copy()
        while len(pf)<8: pf.append("")
        pf[1]=str(new_pid)
        pf[2]=str(new_mid)

        pid_map[old_pid]=new_pid
        add_lines.append("$ Damaged material/property\n")
        add_lines.append(make(mf))
        add_lines.append(make(pf))

    selected_eid_set=set(eid for eid,_ in selected_elems)

    out_lines=[]
    inserted=False

    for line in lines:
        flds=ff(line)
        c=flds[0]

        if not inserted and c=="GRID":
            out_lines.extend(add_lines)
            inserted=True

        if c in elem_cards:
            eid=int(flds[1])
            if eid in selected_eid_set:
                old_pid=int(flds[2])
                if old_pid in pid_map:
                    line=line[:16]+field(pid_map[old_pid])+line[24:]

        out_lines.append(line)

    if not inserted:
        out_lines.extend(add_lines)

    ts=time.strftime("%Y%m%d_%H%M%S")
    out_path=os.path.join(out_dir,f"FEM_{ts}.bdf")

    with open(out_path,"w") as f:
        f.writelines(out_lines)

    print("selected elements:",len(selected_elems))
    print("new properties:",len(pid_map))
    print("pid map:",pid_map)
    print("saved:",out_path)

    return out_path,[eid for eid,_ in selected_elems]

def check_bdf_duplicates(bdf_path):
    from collections import defaultdict

    elem_cards={"CBAR","CBEAM","CROD","CTRIA3","CQUAD4","CTETRA","CHEXA"}
    prop_cards={"PSHELL","PBAR","PBEAM","PROD","PSOLID"}
    mat_cards={"MAT1"}

    ids={"elem":defaultdict(list),"prop":defaultdict(list),"mat":defaultdict(list)}

    with open(bdf_path,"r",errors="ignore") as f:
        for i,line in enumerate(f,1):
            card=line[:8].strip()
            flds=[line[j:j+8].strip() for j in range(0,len(line.rstrip("\n")),8)]

            if card in elem_cards and len(flds)>1 and flds[1].isdigit():
                ids["elem"][int(flds[1])].append(i)
            elif card in prop_cards and len(flds)>1 and flds[1].isdigit():
                ids["prop"][int(flds[1])].append(i)
            elif card in mat_cards and len(flds)>1 and flds[1].isdigit():
                ids["mat"][int(flds[1])].append(i)

    for k,d in ids.items():
        dup={i:lines for i,lines in d.items() if len(lines)>1}
        print(k,"duplicates:",len(dup))
        for i,lines in list(dup.items())[:20]:
            print(i,lines)

    return ids
if __name__ == "__main__":
    import numpy as np
    #有几个node飘在外面
    # print("Parsing BDF file...")
    #
    # nodes, elements = parse_bdf(bdf_path)
    #
    # print(f"Total nodes: {len(nodes)}")
    # print(f"Total elements: {len(elements)}")
    #
    # write_nodes(nodes, nodes_csv)
    # write_elements(elements, elements_csv)
    #
    # print("Done!")
    # print(f"Nodes saved to: {nodes_csv}")
    # print(f"Elements saved to: {elements_csv}")
    # plot_bdf_mesh_3d(data_dir="data")

    # plot_bdf_mesh_3d("data")
    # plot_nodes_fixed_box("data")
    # plot_pick_nodes_fixed_box("data")

    # plot_selected_nodes(1000)
    # check_bdf_duplicates("data/FEM_20260428_202832.bdf")
    selected_nodes = select_nodes_in_ellipsoid(
        node=1000,
        a=30,
        b=30,
        c=30,
        Ome=np.eye(3),
        decay_rate=0.5,
        data_dir="old_data"
    )

    out_bdf, selected_elements = write_damaged_bdf(
        selected_nodes,
        bdf_in="data/FEM_only.bdf",
        damage_ratio=0.1,
        out_dir="old_data"
    )



