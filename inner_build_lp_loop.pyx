from __future__ import division
from cpython cimport array
from cpython cimport mapping
import array
import time
import copy
import sage.all as sg
from sage.all import libgap as lg
import cplex as cp
from libcpp.vector cimport vector
from libcpp.string cimport string
# cimport cplex as ccp

#### Inner LP build loop for cut identifier strings ####
cdef str build_cut_key(list cut_orb):
    cdef str cut_key = ""
    cdef int cut_len = len(cut_orb)
    for i_col in range(cut_len - 1):
        cut_key += str(cut_orb[i_col]) + ","
    cut_key += str(cut_orb[-1])
    return cut_key
#### Inner LP build loop for cut values for full model with cuts added ####
cdef list full_inner_build_vals(list cut_orbs, dict hash_table, dict node_table,
                                dict cut_rhs, double rhs, int num_full_row, 
                                int num_full_cut, int num_nz, dict cut_map):
    cdef vector[int] full_cut_idx
    cdef vector[double] full_cut_val
    cdef int nz_idx = 0
    cdef int cut_len = len(cut_orbs[0])
    cdef int i_col
    cdef double i_val
    cdef int i_orb
    cdef int i_cut
    cut_rows = []
    cut_rows_sense = []
    cut_rows_rhs = []
    cut_rows_name = []
    name = num_full_row
#     print(cut_map)
#     input()
    for i_orb in range(len(cut_orbs)):
#         print(cut_orbs[i_orb])
#         input()
        cut_orb = [cut_map.get(unq) for unq in cut_orbs[i_orb]]
#         print(cut_orb)
#         input()
        in_cut = cut_orb
        in_cut.append(rhs)
        in_cut = [sg.Rational(el) for el in in_cut]
        cut_key = build_cut_key(in_cut)
#         node_key = hash_table.get(cut_key)
        for i_col in range(cut_len):
            i_val = cut_orb[i_col]
            if (abs(i_val - 0) < 1e-6):
                continue
            full_cut_idx.push_back(int(i_col))
            full_cut_val.push_back(float(i_val)) 
            nz_idx += 1 
        cut_rows.append(cp.SparsePair(full_cut_idx, full_cut_val))
        cut_rows_sense.append("G")
        cut_rows_rhs.append(float(rhs))
        cut_rows_name.append("R%d" % name)
        hash_table[cut_key] = num_full_cut
        node_table[num_full_cut] = cut_key
        cut_rhs[num_full_cut] = rhs
        num_full_cut += 1
        full_cut_idx.clear()
        full_cut_val.clear()
    return [cut_rows, cut_rows_sense, cut_rows_rhs, cut_rows_name]
#### Inner LP build loop for cut indices for full model with cuts added ####
cdef vector[int] full_inner_build_indices(list cut_orb, int num_nz):
    cdef vector[int] full_cut_idxs
    cdef int nz_idx = 0
    cdef int cut_len = len(cut_orb)
    cdef int i_col
    cdef double i_val
    for i_col in range(cut_len):
        i_val = cut_orb[i_col]
        if (abs(i_val - 0) < 1e-6):
            continue
        full_cut_idxs.push_back(int(i_col))
        nz_idx += 1 
    return full_cut_idxs
#### Inner ALP build loop for cut values for aggregate model ####
cdef list agg_inner_build_value(dict hash_table, dict node_orbit,
                                 dict col_orbit_rep, int num_col_orbits, int num_row_orbits, 
                                 list first_cut, list cut_orbit_rep, dict col_orbit, 
                                 dict col_orbit_size):
    agg_lin_constrs = []
    cdef vector[double] agg_values
    cdef vector[double] agg_val_sprs
    cdef vector[int] agg_idx_sprs
#     agg_cuts = {agg_row : {} for agg_row in range(num_row_orbits)}
    rhs = first_cut[-1]
    num_col = len(col_orbit)
    cdef int rep_len = len(cut_orbit_rep)
#     for cut_orb in cut_orbit:
#     for i in range(len(cut_orbit)):
#         cut_orb = cut_orbit[i]
#         in_cut = list(cut_orb)
#         in_cut.append(rhs)
#         cut_key = build_cut_key(in_cut)
#         node_key = hash_table.get(cut_key)
#         node_orb = node_orbit.get(node_key)
#         for agg_col in range(num_col_orbits):
#             rep = col_orbit_rep.get(agg_col)
#             agg_cut = agg_cuts.get(node_orb)
#             if agg_cut.get(agg_col) is None:
#                 agg_cuts[node_orb][agg_col] = 0
#             agg_cuts[node_orb][agg_col] += float(cut_orb[rep])
    for i in range(rep_len):
        agg_values.clear()
        agg_values.resize(num_col_orbits)
        cut_orb = cut_orbit_rep[i]
#         in_cut = list(cut_orb)
#         in_cut.append(rhs)
#         cut_key = build_cut_key(in_cut)
#         node_key = hash_table.get(cut_key)
#         node_orb = node_orbit.get(node_key)
        for i_col in range(num_col):
            agg_col = int(col_orbit.get(i_col))
            coeff = cut_orb[i_col]
            if abs(coeff - 0) < 1e-6:
                continue
#             if agg_cuts.get(i).get(agg_col) is None:
#                 agg_cuts[i][agg_col] = 0
#             agg_cuts[i][agg_col] += float(coeff.sage()/col_orbit_size.get(agg_col))
            agg_values[agg_col] += float(coeff/col_orbit_size.get(agg_col))
        for agg_col in range(num_col_orbits):
            val = agg_values[agg_col]
            if (abs(val - 0) < 1e-6):
                continue
            agg_val_sprs.push_back(val)
            agg_idx_sprs.push_back(agg_col)
        agg_lin_constrs.append(cp.SparsePair(agg_idx_sprs, agg_val_sprs))
        agg_idx_sprs.clear()
        agg_val_sprs.clear()
#     print(agg_lin_constrs)
#     input()
    return agg_lin_constrs
#### Inner ALP build loop for cut values for aggregate model single cut####
cdef dict agg_inner_build_value_dict_single(dict hash_table, dict node_orbit, 
                                     dict col_orbit_rep, int num_col_orbits, int num_row_orbits, 
                                     list first_cut, list cut_orbit_rep, dict col_orbit, dict col_orbit_size):
    agg_cuts = {agg_row : {} for agg_row in range(num_row_orbits)}
    rhs = first_cut[-1]
    num_col = len(col_orbit)
    cdef int rep_len = len(cut_orbit_rep)
#     for cut_orb in cut_orbit:
#     for i in range(len(cut_orbit)):
#         cut_orb = cut_orbit[i]
#         in_cut = list(cut_orb)
#         in_cut.append(rhs)
#         cut_key = build_cut_key(in_cut)
#         node_key = hash_table.get(cut_key)
#         node_orb = node_orbit.get(node_key)
#         for agg_col in range(num_col_orbits):
#             rep = col_orbit_rep.get(agg_col)
#             agg_cut = agg_cuts.get(node_orb)
#             if agg_cut.get(agg_col) is None:
#                 agg_cuts[node_orb][agg_col] = 0
#             agg_cuts[node_orb][agg_col] += float(cut_orb[rep])
    for i in range(rep_len):
        cut_orb = cut_orbit_rep[i]
#         in_cut = list(cut_orb)
#         in_cut.append(rhs)
#         cut_key = build_cut_key(in_cut)
#         node_key = hash_table.get(cut_key)
#         node_orb = node_orbit.get(node_key)
        for i_col in range(num_col):
            agg_col = int(col_orbit.get(i_col))
            coeff = cut_orb[i_col]
            if agg_cuts.get(i).get(agg_col) is None:
                agg_cuts[i][agg_col] = 0
            agg_cuts[i][agg_col] += float(coeff/col_orbit_size.get(agg_col))
#     print(agg_cuts)
#     input()
    return agg_cuts
#### Inner ALP build loop for cut rhs for aggregate model
cdef vector[double] agg_inner_build_rhs(dict hash_table, dict node_orbit, 
                                       dict col_orbit_rep, int num_col_orbits, int num_row_orbits, 
                                       list first_cut):
    cdef vector[double] agg_rhs
    agg_rhs.resize(num_row_orbits)
    rhs = first_cut[-1]
    for i in range(num_row_orbits):
#         in_cut = list(cut_orb)
#         in_cut.append(rhs)
#         cut_key = build_cut_key(in_cut)
#         node_key = hash_table.get(cut_key)
#         node_orb = node_orbit.get(node_key)
        agg_rhs[i] = float(rhs)
#     print(agg_rhs)
#     input()
    return agg_rhs

#### Full cython lift cut function ####
def lift_cut_cy_average(cut, group, first_group, orbits, num_tot, hash_table, 
                        node_table, full_mdl, agg_mdl, cut_reps, cut_rhs, agg_cut_start):
    num_full_cut = copy.deepcopy(full_mdl.variables.get_num() + full_mdl.linear_constraints.get_num())
    num_full_row = copy.deepcopy(full_mdl.variables.get_num() + full_mdl.linear_constraints.get_num())
    num_agg_cut = copy.deepcopy(agg_cut_start)
    num_row_orbits = 0
    new_cut_orbit_size = {}
    new_orbits = []
    cut_orbit_rep = []
    new_cut_orbit ={}
    node_orbit = {}
    agg_col_val = {}
    agg_col_idx = {}
    cdef int col
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
    rhs = cut[-1]
    num_nz = 0
    cdef int idx
    cdef int num_unique = 0
    cdef double scale
    cdef int var
    
    cdef int map_val
    unique_coeff = set()
    unique_map = {}
    coeff_map = {}
    for col, coeff in enumerate(cut[:-1]):
        if coeff not in unique_coeff:
            coeff_map[coeff] = num_unique
            unique_map[num_unique] = coeff
            unique_coeff.add(coeff)
            num_unique += 1
        for idx in range(orbits.orbit_col_start[col], orbits.orbit_col_start[col + 1]):
            scale = orbits.col_orbit_size[col]
            var = orbits.orbit_col[idx]
            lifted_cut[var] = coeff
            if abs(coeff - 0) > 1e-6:
                num_nz += 1
    lifted_cut.append(rhs)
    cdef vector[int] cut_for_gap
    for col in range(len(orbits.col_orbit.keys())):
        coeff = lifted_cut[col]
        map_val = coeff_map.get(coeff)
        cut_for_gap.push_back(sg.Integer(map_val))
    lift_str = ",".join(str(el) for el in lifted_cut)
#     print(lift_str)

    # lift_str += "," + str(rhs)
    if hash_table.get(lift_str) != None:
        return (new_orbits, new_cut_orbit_size, [], [], [], [])
    start_1 = time.perf_counter()
    cut_orbit = lg.Orbit(lg(first_group), lg(cut_for_gap), lg.Permuted)
    if len(cut_orbit) > 300000:
        print("too symmetric cuts for this gomory cut...moving on")
        hash_table[lift_str] = "too_big"
        return (new_orbits, new_cut_orbit_size, [], [], [], [])
    for_cut_reps = [unique_map.get(unq) for unq in cut_orbit[0].sage()]
    cut_reps[num_full_cut] = for_cut_reps
#     print(cut_reps)
#     input()
    gen_cut_time = time.perf_counter() - start_1
    print("\n%.2f seconds required to generate %d cuts for full model." % (gen_cut_time, len(cut_orbit)))
    start_2 = time.perf_counter()
    stab_cut_orbit = lg.OrbitsDomain(lg(group), cut_orbit, lg.Permuted)
    orbit_lengths = lg.OrbitLengthsDomain(lg(group), cut_orbit, lg.Permuted)
    num_stab_orbits = len(orbit_lengths)
    gen_cut_agg_orb_time = time.perf_counter() - start_2
    print("%.2f seconds required to generate %d cuts for aggregate model." % (gen_cut_agg_orb_time, len(orbit_lengths)))
    cut_rows = []
    cut_rows_sense = []
    cut_rows_rhs = []
    cut_rows_name = []
    name = num_full_row
    cut_idx = 0
    num_cut = 0
    start_3 = time.perf_counter()
    n_f_cut_added = len(cut_orbit)
    n_a_cut_added = len(orbit_lengths)
    cdef int i_orb
    cdef int i_cut
#     for i_orb in range(len(cut_orbit)):
#         cut_orb = cut_orbit[i_orb]
#         in_cut = cut_orb.sage()
#         in_cut.append(rhs)
#         in_cut = [sg.Rational(el) for el in in_cut]
#         cut_key = build_cut_key(in_cut)
# #         node_key = hash_table.get(cut_key)
#         full_cut_idx = full_inner_build_indices(cut_orb.sage(), num_nz)
#         full_cut_val = full_inner_build_vals(cut_orb.sage(), num_nz)
#         cut_rows.append(cp.SparsePair(full_cut_idx, full_cut_val))
#         cut_rows_sense.append("G")
#         cut_rows_rhs.append(float(rhs))
#         cut_rows_name.append("R%d" % name)
#         hash_table[cut_key] = num_full_cut
#         node_table[num_full_cut] = cut_key
#         cut_rhs[num_full_cut] = rhs
#         num_full_cut += 1
    output = full_inner_build_vals(cut_orbit.sage(), hash_table, node_table,
                                  cut_rhs, rhs, num_full_row, num_full_cut, num_nz,
                                  unique_map)
    cut_rows = output[0]
    cut_rows_sense = output[1]
    cut_rows_rhs = output[2]
    cut_rows_name = output[3]
    full_mdl.linear_constraints.add(lin_expr = cut_rows, senses = cut_rows_sense,
                                       rhs = cut_rows_rhs, names = cut_rows_name)
    add_to_full_model_time += time.perf_counter() - start_3
    print("%.2f seconds required to add cuts to full model." % add_to_full_model_time)
    del cut_orbit
    del cut_rows
    del cut_rows_sense
    del cut_rows_rhs
    del cut_rows_name
    t_5 = time.perf_counter()
    for i_orb in range(num_stab_orbits):
    #     start_2 = time.perf_counter()
        rep_orb = [unique_map.get(unq) for unq in stab_cut_orbit[i_orb][0].sage()]
#         cut_reps[num_full_cut].append(rhs)
        orb_size = orbit_lengths[i_orb]
#         ###### TESTING #######
#         in_cut = list(rep_orb)
#         in_cut.append(rhs)
#         in_cut = [sg.Rational(el) for el in in_cut]
#         cut_key = build_cut_key(in_cut)
#         node_key = hash_table.get(cut_key)
#         if node_orbit.get(node_key) is not None:
#             continue
#         full_cut_idx = full_inner_build_indices(rep_orb.sage(), num_nz)
#         full_cut_val = full_inner_build_vals(rep_orb.sage(), num_nz)
#         cut_rows.append(cp.SparsePair(full_cut_idx, full_cut_val))
#         cut_rows_sense.append("G")
#         cut_rows_rhs.append(float(rhs))
#         cut_rows_name.append("R%d" % name)
#         hash_table[cut_key] = num_full_cut
#         node_table[num_full_cut] = cut_key
#         node_orbit[num_full_cut] = num_row_orbits
#         cut_orbit_rep.append(rep_orb)
#         new_cut_orbit[num_row_orbits] = num_full_cut
#         new_cut_orbit_size[num_agg_cut] = int(orb_size)
#         new_orbits.append([num_full_cut])
#         num_row_orbits += 1
#         num_agg_cut += 1
#         num_full_cut += 1
#         name += 1
        new_orbit = []
        for i_cut in range(orb_size):
#             print(unique_map)
#             print(stab_cut_orbit[i_orb][i_cut])
            cut_orb = [unique_map.get(unq) for unq in stab_cut_orbit[i_orb][i_cut].sage()]
#             print(cut_orb)
#             input()
            in_cut = cut_orb
            in_cut.append(rhs)
            in_cut = [sg.Rational(el) for el in in_cut]
            cut_key = build_cut_key(in_cut)
#             in_cut = list(cut_coeff)
#             in_cut.append(rhs)
#             in_cut = [sg.Rational(el) for el in in_cut]
#             cut_key = build_cut_key(in_cut)
#             node_key = hash_table.get(cut_key)
#             full_cut_idx = full_inner_build_indices(cut_coeff.sage(), num_nz)
#             full_cut_val = full_inner_build_vals(cut_coeff.sage(), num_nz)
#             cut_rows.append(cp.SparsePair(full_cut_idx, full_cut_val))
#             cut_rows_sense.append("G")
#             cut_rows_rhs.append(float(rhs))
#             cut_rows_name.append("R%d" % name)
            cut_node = hash_table.get(cut_key)
#             node_table[num_full_cut] = cut_key
#             cut_rhs[num_full_cut] = rhs
            node_orbit[cut_node] = num_row_orbits 
            new_orbit.append(cut_node)
#             num_full_cut += 1
#             name += 1
        cut_orbit_rep.append(rep_orb)
        new_cut_orbit[num_row_orbits] = new_orbit[0]
        new_cut_orbit_size[num_agg_cut] = int(orb_size)
        new_orbits.append(new_orbit)
        num_row_orbits += 1
        num_agg_cut += 1
#     # LIFT THE CUTS TO THE ORIGINAL MODEL
#     cut_rows = []
#     cut_rows_sense = []
#     cut_rows_rhs = []
#     cut_rows_name = []
#     name = num_full_row
#     cut_idx = 0
#     num_cut = 0
#     for cut_orb in cut_orbit:
#         start_2 = time.perf_counter()
#         ###### TESTING #######
#         in_cut = list(cut_orb)
#         in_cut.append(rhs)
#         in_cut = [sg.RealNumber(el) for el in in_cut]
#         cut_key = build_cut_key(in_cut)
#         node_key = hash_table.get(cut_key)
#         # full_cut_idx = full_inner_build_indices(cut_orb.sage(), num_nz)
#         # full_cut_val = full_inner_build_vals(cut_orb.sage(), num_nz)
#         ##### TESTING #####
#         # cut_rows.append(cp.SparsePair(full_cut_idx, full_cut_val))
#         # cut_rows_sense.append("G")
#         # cut_rows_rhs.append(float(rhs))
#         # cut_rows_name.append("C%d" % name)
#         # name += 1
#         # add_to_full_model_time += time.perf_counter() - start_2
#         if node_orbit.get(node_key) != None:
#             continue
#         new_orbit = []
#         start_3 = time.perf_counter()
#         stab_orb = lg.Orbit(lg(group), lg(cut_orb), lg.Permuted)
# #         print(stab_orb) 
# #         input()
#         for orb in stab_orb:
#             full_cut_idx = full_inner_build_indices(orb.sage(), num_nz)
#             full_cut_val = full_inner_build_vals(orb.sage(), num_nz)
#             cut_rows.append(cp.SparsePair(full_cut_idx, full_cut_val))
#             cut_rows_sense.append("G")
#             cut_rows_rhs.append(float(rhs))
#             cut_rows_name.append("R%d" % name)
#             name += 1
#             in_cut = list(orb)
#             in_cut.append(rhs)
#             in_cut = [sg.RealNumber(el) for el in in_cut]
#             orb_key = build_cut_key(in_cut)
# #             print(orb_key)
# #             input()
#             hash_table[orb_key] = num_full_cut
#             node_table[num_full_cut] = orb_key
#             new_orbit.append(num_full_cut)
#             node_orbit[num_full_cut] = num_row_orbits 
#             num_full_cut += 1
#         cut_orbit_rep.append(stab_orb[0])
#         new_cut_orbit[num_row_orbits] = new_orbit[0]
#         new_cut_orbit_size[num_agg_cut] = len(new_orbit)
#         num_row_orbits += 1
#         num_agg_cut += 1
#         new_orbits.append(new_orbit)
#         gen_cut_agg_orb_time += time.perf_counter() - start_3
# #     full_cut_csr = scipy.sparse.csr_matrix((full_cut_val, full_cut_idx, full_cut_start))
#     full_mdl.linear_constraints.add(lin_expr = cut_rows, senses = cut_rows_sense,
#                                        rhs = cut_rows_rhs, names = cut_rows_name)
    time_to_record_orbits = time.perf_counter() - t_5
    print("%.2f seconds required to record orbits." % time_to_record_orbits)
    del stab_cut_orbit
    # GENERATE AGGREGATE CUTS FOR THE NEW ORBITS IN THE AGGREGATE SPACE
    cut_rows_sense = []
    cut_rows_name = []
    start_4 = time.perf_counter()
    cut_rows = agg_inner_build_value(hash_table, node_orbit, 
                                          orbits.col_orbit_rep, orbits.num_col_orbits, num_row_orbits, 
                                          lifted_cut, cut_orbit_rep, orbits.col_orbit, orbits.col_orbit_size)
    cut_rows_rhs = agg_inner_build_rhs(hash_table, node_orbit, 
                                       orbits.col_orbit_rep, orbits.num_col_orbits, num_row_orbits, 
                                       lifted_cut)
    for i_row in range(num_row_orbits):
        cut_rows_sense.append("G")
#     agg_mdl.linear_constraints.add(lin_expr = cut_rows, senses = cut_rows_sense,
#                                 rhs = cut_rows_rhs, names = cut_rows_name)
    add_to_agg_model_time = time.perf_counter() - start_4
#     print("\n%.2f seconds required to generate %d cuts for full model." % (gen_cut_time, n_f_cut_added))
#     print("%.2f seconds required to generate %d cuts for aggregate model." % (gen_cut_agg_orb_time, num_row_orbits))
    print("%.2f seconds required to add cuts to aggregate model.\n" % add_to_agg_model_time)
#     print(new_orbits, new_cut_orbit_size, cut_rows, cut_rows_sense, cut_rows_rhs, cut_rows_name)
#     input()
    return (new_orbits, new_cut_orbit_size, cut_rows, cut_rows_sense, cut_rows_rhs, cut_rows_name)

#### Full cython lift cut function (only one cut generated) ####
def lift_cut_cy_single(cut, group, first_group, orbits, num_tot, hash_table, 
                        node_table, full_mdl, agg_mdl, cut_reps, cut_rhs, agg_cut_start):
    num_full_cut = copy.deepcopy(num_tot)
    num_full_row = copy.deepcopy(num_tot)
    num_agg_cut = copy.deepcopy(agg_cut_start)
    num_row_orbits = 0
    new_cut_orbit_size = {}
    new_orbits = []
    cut_orbit_rep = []
    new_cut_orbit ={}
    node_orbit = {}
    agg_col_val = {}
    agg_col_idx = {}
    cdef int col
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
    cdef int idx
    cdef int num_unique = 0
    cdef double scale
    cdef int var
    
    cdef int map_val
    unique_coeff = set()
    unique_map = {}
    coeff_map = {}
    for col, coeff in enumerate(cut[:-1]):
        if coeff not in unique_coeff:
            coeff_map[coeff] = num_unique
            unique_map[num_unique] = coeff
            unique_coeff.add(coeff)
            num_unique += 1
        for idx in range(orbits.orbit_col_start[col], orbits.orbit_col_start[col + 1]):
            scale = orbits.col_orbit_size[col]
            var = orbits.orbit_col[idx]
            lifted_cut[var] = coeff
            if abs(coeff - 0) > 1e-6:
                num_nz += 1
    lifted_cut.append(rhs)
    cdef vector[int] cut_for_gap
    for col in range(len(orbits.col_orbit.keys())):
        coeff = lifted_cut[col]
        map_val = coeff_map.get(coeff)
        cut_for_gap.push_back(sg.Integer(map_val))
    lift_str = build_cut_key(lifted_cut)
    # lift_str += "," + str(rhs)
    if hash_table.get(lift_str) != None:
        return (new_orbits, new_cut_orbit_size, [], [], [], [])
    print("Cuts being added: %d" % 1)
    gen_cut_time = time.perf_counter() - start_1
    # LIFT THE CUTS TO THE ORIGINAL MODEL
    cut_rows = []
    cut_rows_sense = []
    cut_rows_rhs = []
    cut_rows_name = []
    name = num_full_row
    start_2 = time.perf_counter()
    ###### TESTING #######
    node_key = hash_table.get(lift_str)
    if node_orbit.get(node_key) != None:
        return new_orbits, new_cut_orbit_size
    new_orbit = []
    start_3 = time.perf_counter()
#     full_cut_idx = full_inner_build_indices(lifted_cut[:-1], num_nz)
#     full_cut_val = full_inner_build_vals(lifted_cut[:-1], num_nz)
#     cut_rows.append(cp.SparsePair(full_cut_idx, full_cut_val))
#     cut_rows_sense.append("G")
#     cut_rows_rhs.append(float(rhs))
#     cut_rows_name.append("R%d" % name)
#     cut_orbit_rep.append(lifted_cut[:-1])
#     name += 1
#     hash_table[lift_str] = num_full_cut
#     node_table[num_full_cut] = lift_str
    cut_reps[num_full_cut] = lifted_cut[:-1]
    cut_rhs[num_full_cut] = rhs
    output = full_inner_build_vals([cut_for_gap], hash_table, node_table,
                                  cut_rhs, rhs, num_full_row, num_full_cut, num_nz,
                                  unique_map)
    cut_rows = output[0]
    cut_rows_sense = output[1]
    cut_rows_rhs = output[2]
    cut_rows_name = output[3]
    new_orbit.append(num_full_cut)
    node_orbit[num_full_cut] = num_row_orbits 
    num_full_cut += 1
    new_cut_orbit[num_row_orbits] = new_orbit[0]
    new_cut_orbit_size[num_agg_cut] = len(new_orbit)
    cut_orbit_rep.append(lifted_cut[:-1])
    num_row_orbits += 1
    num_agg_cut += 1
    new_orbits.append(new_orbit)
    gen_cut_agg_orb_time += time.perf_counter() - start_3
    start_2 = time.perf_counter()
    full_mdl.linear_constraints.add(lin_expr = cut_rows, senses = cut_rows_sense,
                                       rhs = cut_rows_rhs, names = cut_rows_name)
    add_to_full_model_time += time.perf_counter() - start_2
    # GENERATE AGGREGATE CUTS FOR THE NEW ORBITS IN THE AGGREGATE SPACE
    start_4 = time.perf_counter()
    cut_rows = []
    cut_rows_sense = []
    cut_rows_name = []
    """
    agg_cuts = {agg_row : {} for agg_row in range(num_row_orbits)}
    agg_cuts_rhs = {agg_row : 0 for agg_row in range(num_row_orbits)}
    for cut_orb in cut_orbit:
        in_cut = list(cut_orb)
        in_cut.append(rhs)
        cut_key = build_cut_key(in_cut)
        node_key = hash_table.get(cut_key)
        node_orb = node_orbit.get(node_key)
        for agg_col in range(orbits.num_col_orbits):
            rep = orbits.col_orbit_rep.get(agg_col)
            if agg_cuts.get(node_orb).get(agg_col) is None:
                agg_cuts[node_orb][agg_col] = 0
            agg_cuts[node_orb][agg_col] += float(cut_orb[rep])
        agg_cuts_rhs[node_orb] += float(rhs)
    print(agg_cuts)
    print(agg_cuts_rhs)
    input()
    """
    cut_rows = agg_inner_build_value(hash_table, node_orbit, 
                                          orbits.col_orbit_rep, orbits.num_col_orbits, num_row_orbits, 
                                          lifted_cut, cut_orbit_rep, orbits.col_orbit, orbits.col_orbit_size)
    cut_rows_rhs = agg_inner_build_rhs(hash_table, node_orbit, 
                                       orbits.col_orbit_rep, orbits.num_col_orbits, num_row_orbits, 
                                       lifted_cut)
#     for agg_cut, vec in agg_cuts.items():
#         agg_cut_idx = []
#         agg_cut_val = []
#         for agg_col, val in vec.items():
#             agg_cut_idx.append(agg_col)
#             agg_cut_val.append(val)
#         cut_rows.append(cp.SparsePair(agg_cut_idx, agg_cut_val))
#         cut_rows_sense.append("G")
#         cut_rows_rhs.append(agg_cuts_rhs.get(agg_cut)) 
        #cut_rows_name.append("C%d" % node_orb)
    for i_row in range(num_row_orbits):
        cut_rows_sense.append("G")
#     agg_mdl.linear_constraints.add(lin_expr = cut_rows, senses = cut_rows_sense,
#                                 rhs = cut_rows_rhs, names = cut_rows_name)
    add_to_agg_model_time = time.perf_counter() - start_4
    print("\n%.2f seconds required to generate %d cuts for full model." % (gen_cut_time, 1))
    print("%.2f seconds required to generate %d cuts for aggregate model." % (gen_cut_agg_orb_time, num_row_orbits))
    print("%.2f seconds required to add cuts to full model." % add_to_full_model_time)
    print("%.2f seconds required to add cuts to aggregate model.\n" % add_to_agg_model_time)
    return (new_orbits, new_cut_orbit_size, cut_rows, cut_rows_sense, cut_rows_rhs, cut_rows_name)

#### Aggregate the A matrix to get the sum aggregate formulation #####
def aggregate_A_mat_cy(mdl, orbits):
    agg_mdl = cp.Cplex()
    num_col_orbits, num_row_orbits = orbits.get_num_orbits()
    col_orbit, row_orbit = orbits.get_orbits()
    col_orbit_rep, row_orbit_rep = orbits.get_orbit_reps()
    col_orbit_size, row_orbit_size = orbits.get_orbit_sizes()
    # ADD AGGREGATE ROWS
    rhs = mdl.linear_constraints.get_rhs()
    names = mdl.linear_constraints.get_names()
    sense = mdl.linear_constraints.get_senses()
    agg_rhs = []
    agg_name = []
    agg_sense = []
    for agg_row in range(num_row_orbits):
        rep = row_orbit_rep.get(agg_row)
        agg_rhs.append(rhs[rep] * row_orbit_size[agg_row])
        agg_sense.append(sense[rep])
        agg_name.append(names[rep])
    agg_mdl.linear_constraints.add(rhs = agg_rhs, senses = agg_sense, names = agg_name)
    # BUILD AGGREGATE COLUMNS
    objs = mdl.objective.get_linear()
    lbs = mdl.variables.get_lower_bounds()
    ubs = mdl.variables.get_upper_bounds()
#     types = mdl.variables.get_types()
    names = mdl.variables.get_names()
    agg_obj = []
    agg_col_vec = []
    agg_lb = []
    agg_ub = []
    agg_type = []
    agg_name = []
    for agg_col in range(num_col_orbits):
        agg_value = {}
        rep = col_orbit_rep.get(agg_col)
        index, value = mdl.variables.get_cols(rep).unpack()
        for ind, val in zip(index, value):
            row_orb = row_orbit.get(ind)
            if agg_value.get(row_orb) == None:
                agg_value[row_orb] = 0
            agg_value[row_orb] += val
        row_val = []
        row_idx = []
#         print(agg_value)
#         input()
        for agg_row, val in agg_value.items():
            row_val.append(val)
            row_idx.append(int(agg_row))
        agg_obj.append(objs[rep])
        agg_col_vec.append(cp.SparsePair(row_idx, row_val))
        agg_lb.append(lbs[rep] * col_orbit_size[agg_col])
        agg_ub.append(ubs[rep] * col_orbit_size[agg_col])
        agg_type.append("C")
        agg_name.append(names[rep])
    agg_mdl.variables.add(obj = agg_obj, lb = agg_lb, ub = agg_ub, types = agg_type,
                        names = agg_name, columns = agg_col_vec)
#     agg_mdl.write("agg_model_cplex.lp")
    return agg_mdl
#### Aggregate the A mat to get the average aggregate formulation ####
def aggregate_A_mat_average_cy(mdl, orbits, timer):
    cdef int time_out = 0
    t_0 = time.perf_counter()
    t_1 = time.perf_counter()
    agg_mdl = cp.Cplex()
    cdef int num_col_orbits
    cdef int num_row_orbits
    num_col_orbits, num_row_orbits = orbits.get_num_orbits()
    col_orbit, row_orbit = orbits.get_orbits()
    col_orbit_rep, row_orbit_rep = orbits.get_orbit_reps()
    col_orbit_size, row_orbit_size = orbits.get_orbit_sizes()
    num_cols = len(orbits.col_orbit)
    # BUILD AGGREGATE COLUMNS
#     objs = mdl.objective.get_linear()
#     lbs = mdl.variables.get_lower_bounds()
#     ubs = mdl.variables.get_upper_bounds()
#     names = mdl.variables.get_names()
    cdef vector[double] agg_obj
    cdef vector[double] agg_lb
    cdef vector[double] agg_ub
    cdef int agg_row
    cdef int i_mat
    cdef int agg_col
    cdef int i_col
    cdef int rep
    cdef double bound_scale
    cdef double i_val
#     cdef vector[string] agg_type
#     cdef vector[string] agg_name
#     agg_obj = []
#     agg_lb = []
#     agg_ub = []
    agg_type = []
    agg_name = []
    for agg_col in range(num_col_orbits):
        bnd_scale = col_orbit_size[agg_col]
        rep = col_orbit_rep.get(agg_col)
        agg_obj.push_back(mdl.objective.get_linear(rep))
        agg_lb.push_back(mdl.variables.get_lower_bounds(rep) * bnd_scale)
        agg_ub.push_back(mdl.variables.get_upper_bounds(rep) * bnd_scale)
        agg_type.append("C")
        agg_name.append(mdl.variables.get_names(rep))
    agg_mdl.variables.add(obj = agg_obj, lb = agg_lb, ub = agg_ub, types = agg_type,
                          names = agg_name)
    print("Time to add aggregate variables: %.2f" % (time.perf_counter() - t_1))
    t_1 = time.perf_counter()
    # ADD AGGREGATE ROWS
#     rhs = mdl.linear_constraints.get_rhs()
#     names = mdl.linear_constraints.get_names()
#     sense = mdl.linear_constraints.get_senses()
    cdef vector[double] agg_rhs
#     agg_rhs = []
    agg_name = []
    agg_sense = []
    agg_lin_constrs = []
    cdef vector[int] idx
    cdef vector[double] val
    cdef vector[double] agg_values
    cdef vector[double] agg_val_sprs
    cdef vector[int] agg_idx_sprs
    cdef int num_iter
    for agg_row in range(num_row_orbits):
        if timer[0].get_run_time() > 7200:
            time_out = 1
            break
        rep = row_orbit_rep.get(agg_row)
        agg_rhs.push_back(mdl.linear_constraints.get_rhs(rep))
        agg_sense.append(mdl.linear_constraints.get_senses(rep))
        agg_name.append(mdl.linear_constraints.get_names(rep))
        idx, val = mdl.linear_constraints.get_rows(rep).unpack()
        num_iter = len(idx)
#         cdef vector[int] agg_indices
#         agg_indices.resize(num_col_orbits)
        agg_values.clear()
        agg_values.resize(num_col_orbits)
#         agg_row = {}
        for i_mat in range(num_iter):
            i_col = idx[i_mat]
            i_val = val[i_mat]
            agg_col = int(orbits.col_orbit.get(i_col))
            scale = col_orbit_size[agg_col]
#             if agg_row.get(agg_col) is None:
#                 agg_row[agg_col] = 0
#             agg_row[agg_col] += i_val/scale
            agg_values[agg_col] += i_val/scale
        for agg_col in range(num_col_orbits):
            i_val = agg_values[agg_col]
            if (abs(i_val - 0) < 1e-6):
                continue
            agg_val_sprs.push_back(i_val)
            agg_idx_sprs.push_back(agg_col)
#         agg_idx = list(key for key in agg_row.keys())
#         agg_val = list(agg_row.get(key) for key in agg_row.keys())
        agg_lin_constrs.append(cp.SparsePair(agg_idx_sprs, agg_val_sprs))
        agg_idx_sprs.clear()
        agg_val_sprs.clear()
    agg_mdl.linear_constraints.add(lin_expr = agg_lin_constrs, rhs = agg_rhs, senses = agg_sense, names = agg_name)
#     agg_mdl.write("agg_model_cplex.lp")
#     input()
    build_time = time.perf_counter() - t_0
    print("Time to add aggregate constraints: %.2f" % (time.perf_counter() - t_1))
    print("Time to build aggregate model before cuts: %.2f" % build_time)
    return agg_mdl, time_out

def append_level_swapping_cuts(cut_orbits, agg_mdl, mdl, orbits):
    num_col_orbits, num_row_orbits = orbits.get_num_orbits()
    col_orbit, row_orbit = orbits.get_orbits()
    col_orbit_rep, row_orbit_rep = orbits.get_orbit_reps()
    col_orbit_size, row_orbit_size = orbits.get_orbit_sizes()
    cdef int i_orb
    cdef int rep
    cdef int agg_col
    cdef int i_mat
    cdef int i_col
    cdef double i_val
    cdef double scale
    rhs = mdl.linear_constraints.get_rhs()
    names = mdl.linear_constraints.get_names()
    sense = mdl.linear_constraints.get_senses()
    cdef int num_col = mdl.variables.get_num()
    cdef vector[double] agg_rhs
    agg_name = []
    agg_sense = []
    agg_lin_constrs = []
    cdef vector[double] agg_values
    cdef vector[double] agg_val_sprs
    cdef vector[int] agg_idx_sprs
    for i_orb in range(len(cut_orbits)):
        rep = cut_orbits[i_orb][0] - num_col
        agg_rhs.push_back(rhs[rep])
        agg_sense.append(sense[rep])
        agg_name.append(names[rep])
        idx, val = mdl.linear_constraints.get_rows(rep).unpack()
        num_iter = len(idx)
        agg_values.clear()
        agg_values.resize(num_col_orbits)
#         agg_row = {}
        for i_mat in range(num_iter):
            i_col = idx[i_mat]
            i_val = val[i_mat]
            agg_col = int(orbits.col_orbit.get(i_col))
            scale = col_orbit_size[agg_col]
#             if agg_row.get(agg_col) is None:
#                 agg_row[agg_col] = 0
#             agg_row[agg_col] += i_val/scale
            agg_values[agg_col] += i_val/scale
        for agg_col in range(num_col_orbits):
            i_val = agg_values[agg_col]
            if (abs(i_val - 0) < 1e-6):
                continue
            agg_val_sprs.push_back(i_val)
            agg_idx_sprs.push_back(agg_col)
#         agg_idx = list(key for key in agg_row.keys())
#         agg_val = list(agg_row.get(key) for key in agg_row.keys())
#         print(agg_idx_sprs)
#         print(agg_val_sprs)
#         input()
        agg_lin_constrs.append(cp.SparsePair(agg_idx_sprs, agg_val_sprs))
        agg_idx_sprs.clear()
        agg_val_sprs.clear()
        
        
    
    