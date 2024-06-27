
import MDAnalysis as mda
import numpy as np







def get_vectors(lista
                ):
    r"""This function gets a list with a specific carbon (e.g. C34 or C22) 
    and :program: its respective hidrogens (e.g. H4X, H4Y). It computes the vectors
    that connect the carbons and the hydrogens, computes the :math:`cos(\theta)^2` 
    of the angle of these vectors and return the mean of it :math:`\braket{cos(\theta)^2}`

    Parameters
    ----------
    lista : list 
        Vector of the shape :math:`[C*i, HiX,HiY, HiZ]`, the minimun len is 2 (when the carbon
        only have one hydrogen) and the maximun is 4 (when there is three hydrogens) 
        Note: If there is N lipids, there will be N carbons and the i represents
        the position of the carbon in the lipid tail

    Returns
    -------
    order : float 
        Float with the mean of :math:`\braket{cos(\theta)^2}`

    Notes
    -----
    The average of the angle of the i-th carbon for all the lipids in the selection is computed
    as follows:

    .. math:: \braket{cos(\theta_i)^2}

    where :math:`\theta_i` is the angle between the vector that connect the i-th carbon and the hydrogen.

    

    """
    angles = []
    for i in range(len(lista)-1):
        vectores =lista[i+1].positions - lista[0].positions
        angles.append(vectores)
    angles = np.concatenate(angles, axis = 0)
    costheta = angles[:,2]**2/np.linalg.norm(angles, axis=1)**2
    order = np.mean(costheta)
    return order

def order_sn1(sel, lipid, n_chain):

    r"""

    Code to loop over the number of carbons in the lipid tail and return a vector 
    :math:`[C3i, HiX, HiY, ...]` that is used in get_vectors
    Then, make an array of dim n_chain with the mean :math:`\braket{cos(\theta_i)^2}`

    Parameters
    ----------
    lipid : str 
        Name of the lipid to compute the order parameters
    n_chain : int 
        Number of carbond in the lipid tail sn1

    Returns
    -------
    chains : ndarray 
        Vector dim n_chains with the mean 
    
    Notes
    -----   
    The return is a vector containing the value :math:`\braket{cos^2(\theta_i}`. As follows:

    .. math:: [\braket{cos^2(\theta_2}, \braket{cos^2(\theta_3}, ..., \braket{cos^2(\theta_{n_chain}}]

    The index starts at 2 because since that carbond we account for the lipid tail. 

    """

    chains = []
    for i in range(n_chain):
        selections = [
                f"name C3{i+2}",
                f"name H{i+2}X and not name HX",
                f"name H{i+2}Y and not name HY",
                f"name H{i+2}Z and not name HZ",
        ]
        lista = []
        for selection in selections:
            atoms = sel.select_atoms(selection)
            if atoms.n_atoms != 0:
                lista.append(atoms)
        chains.append(get_vectors(lista))
    chains = np.array(chains)

    return chains


def order_sn2(sel, lipid, n_chain):
    
    r"""

    Code to loop over the number of carbons in the lipid tail and return a vector 
    :math:`[C3i, HiX, HiY, ...]` that is used in get_vectors
    Then, make an array of dim n_chain with the mean :math:`\braket{cos(\theta_i)^2}`

    Parameters
    ----------
    lipid : str 
        Name of the lipid to compute the order parameters
    n_chain : int 
        Number of carbond in the lipid tail sn1

    Returns
    -------
    chains : ndarray 
        Vector dim n_chains with the mean 
    
    Notes
    -----   
    The return is a vector containing the value :math:`\braket{cos^2(\theta_i}`. As follows:

    .. math:: [\braket{cos^2(\theta_2}, \braket{cos^2(\theta_3}, ..., \braket{cos^2(\theta_{n_chain}}] 

    The index starts at 2 because since that carbond we account for the lipid tail. 
    """
    chains = []

    for i in range(n_chain):
        selections = [
                f"name C2{i+2}",
                f"name H{i+2}R and not name HR",
                f"name H{i+2}S and not name HS",
                f"name H{i+2}T and not name HT",
        ]
        if lipid == "POPE" or lipid == "POPS":
            if selections[0] == "name C29":
                selections[1] == "name H91"
            if selections[0] == "name C210":
                selections[1] == "name H101"
        lista = []
        for selection in selections:
            atoms = sel.select_atoms(selection)
            if atoms.n_atoms != 0:
                lista.append(atoms)
        chains.append(get_vectors(lista))
    chains = np.array(chains)
    return chains


def sn1(u, sel_string, lipid, n_chain, start = 0, stop = -1, step = 1):
    r"""
    
    Code to compute the order parameters for the sn1 tail. 
    This sn1 tail is caracterized by C3* as the name of the carbons and H*X, H*Y, H*Z for the hydrogens
    

    Parameters
    ----------
    u : universe 
        universe element or selection element
    sel_string : str 
        selection of lipids to compute SCD (e.g "resname DSPC")
    lipid : str
        lipid to compute the SCD. (May be modified in the future since is a bit redundant)
    n_chain : int 
        integer containing the number of carbons in the lipid tail
    start : int 
        initial frame to start to compute the SCD
    stop : int 
        final frame to compute SCD
    step : int
        skip frames



    Returns
    -------
    order_parameters : (array dim n_chain) 
        array with the values of SCD
    serror : (array dim n_chain) 
        array with the std error of the values
    

    Notes
    -----
    The formula to compute order parameters is

    .. math:: \frac{1}{2} |\braket{3cos(\theta)^2-1}|
    """
    order_parameters = []
    for ts in u.trajectory[start:stop:step]:
        sel = u.select_atoms(sel_string)
        order_parameters.append(order_sn1(sel, lipid, n_chain))
    order_parameters = np.array(order_parameters)
    serror = 1.5 * np.std(order_parameters, axis = 0)/np.sqrt(len(order_parameters))
    order_parameters = np.abs(1.5 * np.mean(order_parameters, axis = 0) - 1)
    return order_parameters, serror



def sn2(u, sel_string, lipid, n_chain, start = 0, stop = -1, step = 1):
    r"""
    
    Code to compute the order parameters for the sn2 tail. 
    This sn1 tail is caracterized by C2* as the name of the carbons and H*X, H*Y, H*Z for the hydrogens
    

    Parameters
    ----------
    u : universe 
        universe element or selection element
    sel_string : str 
        selection of lipids to compute SCD (e.g "resname DSPC")
    lipid : str
        lipid to compute the SCD. (May be modified in the future since is a bit redundant)
    n_chain : int 
        integer containing the number of carbons in the lipid tail
    start : int 
        initial frame to start to compute the SCD
    stop : int 
        final frame to compute SCD
    step : int
        skip frames



    Returns
    -------
    order_parameters : (array dim n_chain) 
        array with the values of SCD
    serror : (array dim n_chain) 
        array with the std error of the values
    

    Notes
    -----
    The formula to compute order parameters is

    .. math:: \frac{1}{2} |\braket{3cos(\theta)^2-1}|
    """  
    order_parameters = []
    for ts in u.trajectory[start:stop:step]:
        sel = u.select_atoms(sel_string)
        order_parameters.append(order_sn2(sel, lipid, n_chain))
    order_parameters = np.array(order_parameters)
    serror = 1.5 * np.std(order_parameters, axis = 0)/np.sqrt(len(order_parameters))
    order_parameters = np.abs(1.5 * np.mean(order_parameters, axis = 0) - 1)
    return order_parameters, serror


