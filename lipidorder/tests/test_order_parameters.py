from lipidorder import lipid_order
import MDAnalysis as mda
import matplotlib.pyplot as plt
import numpy as np
import pytest
import os


@pytest.fixture
def universe():
    gro = "membrane.gro"
    xtc = "membrane.xtc"
    path = os.getcwd()
    last_directory = (os.path.basename(path))
    if last_directory == "lipidorderkit":
        path = os.path.join(path,"lipidorder", "tests")
    elif last_directory == "lipidorder":
        path = os.path.join(path, "tests")
    
    return mda.Universe(os.path.join(path, gro), os.path.join(path,xtc))

@pytest.fixture
def n_chain():
    return 17

@pytest.fixture
def data_test():
    data = {
            'start': [i * 5 for i in range(5) ],
            'stop': [i * 10 for i in range(5)],
            'lipids' : ["DSPC", "POPS"],
        }

    lipids = {
            "DSPC" : [17,17],
            "POPS":[15,17],
        }

    return data, lipids

def test_len_lipidorder(universe, n_chain):
    """Sample test, will always pass so long as import statement worked"""
    print(universe)
    u = universe
    sn1 = lipid_order.sn1(u, "resname POPS", "POPS", 15)
    sn2 = lipid_order.sn2(u, "resname POPS", "POPS",17)
    assert len(sn1) == 2 and len(sn2) == 2

def test_get_angles(universe):
    u = universe
    lipid = "POPS"

    lipids_sel = u.select_atoms(f"resname {lipid}")
    print(lipids_sel)
    carbons = [lipids_sel.select_atoms("name C29"),
               lipids_sel.select_atoms("name H91"),

               ]
    print(carbons)
    angle = lipid_order.get_vectors(carbons)
    assert type(angle) == np.float32

def test_lenchain_lipidorder(universe, n_chain):
    """Sample test, will always pass so long as import statement worked"""
    u = universe
    sn1 = lipid_order.sn1(u, "resname DSPC", "DSPC", n_chain)
    sn2 = lipid_order.sn2(u, "resname DSPC", "DSPC", n_chain)
    assert len(sn1[0]) == n_chain and len(sn2[0]) == n_chain


def test_lipids(universe, data_test):
    data, lipids = data_test
    u = universe
    for start in data["start"]:
        for stop in data["stop"]:
            if stop >start:
                for lipid in data["lipids"]:
                        sn1 = lipid_order.sn1(u, f"resname {lipid}", lipid, lipids[lipid][0])
                        sn2 = lipid_order.sn2(u, f"resname {lipid}", lipid, lipids[lipid][1])
                        assert len(sn1[0]) == lipids[lipid][0] and len(sn2[0]) == lipids[lipid][1] 









"""
n = len(u.trajectory[500:])
block_size = [5*i for i in range(1,15)]
error = []
means = []
for size in block_size:
    blocks = n // size
    print(blocks, size)
    averages = []
    for i in range(blocks):
        print(i*size , (i+1)*size)
        sn1 = lipid_order.sn1(u, "resname DSPC", "DSPC", 17, start = 500+i * size, stop =500+ (i+1) * size)
        averages.append(sn1[0])
    averages = np.array(averages)
    mean = np.mean(averages, axis = 0)
    std = np.var(averages,ddof = 1, axis = 0)
    std = np.sqrt(std/blocks)
    print(f"shape std {std.shape}")
    print(f"shape mean {mean.shape}")

    error.append(std)
    means.append(mean)


means = np.array(means)
error = np.array(error) # expected len(block_size[])X17
c_0 = error[:,0]

print(f"shape std {error.shape}")
print(f"shape mean {means.shape}")


print(n)

print(sn1)

plt.errorbar(block_size, means[:,0] ,yerr = c_0 )
plt.show()
"""



