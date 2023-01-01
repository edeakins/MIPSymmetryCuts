from cpython cimport array
import array
import time
import copy
import sage.all as sg
from sage.all import libgap as lg
import cplex as cp

#### Inner LP build loop for cut identifier strings ####
cdef str build_cut_key(list cut_orb):
    cdef str cut_key = ""
    cdef int cut_len = len(cut_orb)
    for i_col in range(cut_len - 1):
        cut_key += str(cut_orb[i_col]) + ","
    cut_key += str(cut_orb[-1])
    return cut_key
#### Inner LP build loop for cut values for full model with cuts added ####
cdef list full_inner_build_vals(list cut_orb, int num_nz):
    cdef array.array full_cut_vals = array.array('d', [0 for i in range(num_nz)])
    cdef double[:] full_cut_val = full_cut_vals
    cdef int nz_idx = 0
    cdef int cut_len = len(cut_orb)
    cdef int i_col
    cdef double i_val
    for i_col in range(cut_len):
        i_val = cut_orb[i_col]
        if (abs(i_val - 0) < 1e-6):
            continue
        full_cut_val[nz_idx] = float(i_val) 
        nz_idx += 1 
    return full_cut_vals.tolist()
#### Inner LP build loop for cut indices for full model with cuts added ####
cdef list full_inner_build_indices(list cut_orb, int num_nz):
    cdef array.array full_cut_idxs = array.array('i', [-1 for i in range(num_nz)])
    cdef int[:] full_cut_idx = full_cut_idxs
    cdef int nz_idx = 0
    cdef int cut_len = len(cut_orb)
    cdef int i_col
    cdef double i_val
    for i_col in range(cut_len):
        i_val = cut_orb[i_col]
        if (abs(i_val - 0) < 1e-6):
            continue
        full_cut_idx[nz_idx] = int(i_col) 
        nz_idx += 1 
    return full_cut_idxs.tolist()
#### Full cython lift cut function ####
def lift_cut_cy(cut, group, first_group, orbits, num_tot, hash_table, node_table, full_mdl, agg_mdl):
    num_full_cut = copy.deepcopy(num_tot)
    num_full_row = copy.deepcopy(num_tot)
    num_row_orbits = 0
    new_cut_orbit_size = {}
    new_orbits = []
    new_cut_orbit ={}
    node_orbit = {}
    agg_col_val = {}
    agg_col_idx = {}
    lifted_cut = [0 for col in orbits.col_orbit.keys()]
    # TESTING TIMERS
    gen_cut_time = 0
    start_1 = 0
    gen_cut_agg_orb_time = 0
    start_2 = 0
    add_to_full_model_time = 0
    start_3 = 0
    add_to_agg_model_time = 0
    start_4 = 0
    # GENERATE ALL CUTS IN THE ORBIT OF ORIGINAL CUT USING ORIGINAL GROUP G
    start_1 = time.perf_counter()
    rhs = cut[-1]
    num_nz = 0
    for col, coeff in enumerate(cut[:-1]):
        for idx in range(orbits.orbit_col_start[col], orbits.orbit_col_start[col + 1]):
            var = orbits.orbit_col[idx]
            lifted_cut[var] = coeff
            if abs(coeff - 0) > 1e-6:
                num_nz += 1
    lifted_cut.append(rhs)
    lift_str = ",".join(str(el) for el in lifted_cut)
    # lift_str += "," + str(rhs)
    if hash_table.get(lift_str) != None:
        return (None, None, new_orbits, new_cut_orbit_size)
    cut_orbit = lg.Orbit(lg(first_group), lifted_cut[:-1], lg.Permuted)
    gen_cut_time = time.perf_counter() - start_1
    # LIFT THE CUTS TO THE ORIGINAL MODEL
    cut_rows = []
    cut_rows_sense = []
    cut_rows_rhs = []
    cut_rows_name = []
    name = num_full_row
    for cut_orb in cut_orbit:
        start_2 = time.perf_counter()
        ###### TESTING #######
        in_cut = list(cut_orb)
        in_cut.append(rhs)
        cut_key = build_cut_key(in_cut)
        node_key = hash_table.get(cut_key)
        full_cut_idx = full_inner_build_indices(cut_orb.sage(), num_nz)
        full_cut_val = full_inner_build_vals(cut_orb.sage(), num_nz)
        ##### TESTING #####
        cut_rows.append(cp.SparsePair(full_cut_idx, full_cut_val))
        cut_rows_sense.append("G")
        cut_rows_rhs.append(float(rhs))
        cut_rows_name.append("C%d" % name)
        name += 1
        add_to_full_model_time += time.perf_counter() - start_2
        if node_orbit.get(node_key) != None:
            continue
        new_orbit = []
        start_3 = time.perf_counter()
        stab_orb = lg.Orbit(lg(group), cut_orb, lg.Permuted)
        for orb in stab_orb:
            orb_key = ",".join(str(el) for el in orb)
            orb_key += "," + str(lg(rhs))
            hash_table[orb_key] = num_full_cut
            node_table[num_full_cut] = orb_key
            new_orbit.append(num_full_cut)
            node_orbit[num_full_cut] = num_row_orbits 
            num_full_cut += 1
        new_cut_orbit[num_row_orbits] = new_orbit[0]
        new_cut_orbit_size[num_row_orbits] = len(new_orbit)
        num_row_orbits += 1
        new_orbits.append(new_orbit)
        gen_cut_agg_orb_time += time.perf_counter() - start_3
#     full_cut_csr = scipy.sparse.csr_matrix((full_cut_val, full_cut_idx, full_cut_start))
    start_2 = time.perf_counter()
    full_mdl.linear_constraints.add(lin_expr = cut_rows, senses = cut_rows_sense,
                                       rhs = cut_rows_rhs, names = cut_rows_name)
    add_to_full_model_time += time.perf_counter() - start_2
    # GENERATE AGGREGATE CUTS FOR THE NEW ORBITS IN THE AGGREGATE SPACE
    start_4 = time.perf_counter()
    cut_rows = []
    cut_rows_sense = []
    cut_rows_rhs = []
    cut_rows_name = []
    agg_cuts = {agg_row : {} for agg_row in range(num_row_orbits)}
    agg_cuts_rhs = {agg_row : 0 for agg_row in range(num_row_orbits)}
    for cut_orb in cut_orbit:
        cut_key = ",".join(str(el) for el in cut_orb)
        cut_key += "," + str(rhs)
        node_key = hash_table.get(cut_key)
        node_orb = node_orbit.get(node_key)
        for agg_col in range(orbits.num_col_orbits):
            rep = orbits.col_orbit_rep.get(agg_col)
            if agg_cuts.get(node_orb).get(agg_col) is None:
                agg_cuts[node_orb][agg_col] = 0
            agg_cuts[node_orb][agg_col] += float(cut_orb[rep])
        agg_cuts_rhs[node_orb] += float(rhs)
    for agg_cut, vec in agg_cuts.items():
        agg_cut_idx = []
        agg_cut_val = []
        for agg_col, val in vec.items():
            agg_cut_idx.append(agg_col)
            agg_cut_val.append(val)
        cut_rows.append(cp.SparsePair(agg_cut_idx, agg_cut_val))
        cut_rows_sense.append("G")
        cut_rows_rhs.append(agg_cuts_rhs.get(agg_cut)) 
        cut_rows_name.append("C%d" % node_orb)
    agg_mdl.linear_constraints.add(lin_expr = cut_rows, senses = cut_rows_sense,
                                rhs = cut_rows_rhs, names = cut_rows_name)
    add_to_agg_model_time = time.perf_counter() - start_4
    print("\n%.2f seconds required to generate %d cuts for full model." % (gen_cut_time, len(cut_orbit)))
    print("%.2f seconds required to generate %d cuts for aggregate model." % (gen_cut_agg_orb_time, num_row_orbits))
    print("%.2f seconds required to add cuts to full model." % add_to_full_model_time)
    print("%.2f seconds required to add cuts to aggregate model.\n" % add_to_agg_model_time)
    input()
    return (new_orbits, new_cut_orbit_size)