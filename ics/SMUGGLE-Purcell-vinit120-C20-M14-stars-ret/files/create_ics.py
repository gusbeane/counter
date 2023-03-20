import arepo
import numpy as np
import os
import sys

from numba import njit

R_Sgr = 80.
V_Sgr = 120.
pro = -1.

Sgr_center = np.array([R_Sgr, 0., 0.])
Sgr_vel = np.array([0., pro*V_Sgr, 0.])

Rcut = 4.0
BoxSize = 400.0
center = np.array([BoxSize/2.0, BoxSize/2.0, BoxSize/2.0])

Nbody_snapnum = 300

def Menc(Rcut, sn):
    ans = 0.0
    for i in range(6):
        if sn.NumPart_Total[i] > 0:
            part = getattr(sn, 'part'+str(i))
            r = np.linalg.norm(part.pos, axis=1)
            key = np.where(r < Rcut)[0]
            if sn.MassTable[i] > 0.0:
                ans += sn.MassTable[i] * len(key)
            else:
                ans += part.mass[key]

    return ans

def rt_eq2(R_Sgr, sn_MW, sn_Sgr):
    M_Sgr = 0.0
    M_MW = 0.0

    for i in range(6):
        M_MW += sn_MW.MassTable[i] * sn_MW.NumPart_Total[i]
        M_Sgr += sn_Sgr.MassTable[i] * sn_Sgr.NumPart_Total[i]

    rt = R_Sgr * (M_Sgr / (3*M_MW))**(1./3.)

    return rt

def rt_eq6(R_Sgr, sn_MW, sn_Sgr):
    M_MW = Menc(R_Sgr, sn_MW)

    # for a hernquist halo (the only material at R_Sgr)
    rt_guess = rt_eq2(R_Sgr, sn_MW, sn_Sgr)

    for _ in range(10):
        rt = R_Sgr * ((Menc(rt_guess, sn_Sgr)/M_MW))**(1./3.)
        rt_guess = rt

    return rt

def compute_tidal_radius(R_Sgr, sn_MW, sn_Sgr):
    return rt_eq6(R_Sgr, sn_MW, sn_Sgr)

sn = arepo.Snapshot('./', Nbody_snapnum, combineFiles=True)
sn_gas = arepo.Snapshot('MW_ICs.dat-with-grid.hdf5')
sn_Sgr = arepo.Snapshot('Sgr_ICs.dat')

# create extra part type for Sgr stars
npart = np.array( [0,  0,  0,  0,  0,  0,  0])
masses = [0., 0., 0., 0., 0., 0., 0.]
for i in range(6):
    npart[i] = sn.NumPart_Total[i]
    masses[i] = sn.MassTable[i]

npart[0] = sn_gas.NumPart_Total[0]
masses[0] = sn_gas.MassTable[0]

npart[5] = sn_Sgr.NumPart_Total[1]
masses[5] = sn_Sgr.MassTable[1]

# now create vacuum
pos_gas = sn_gas.part0.pos - center
R = np.linalg.norm(pos_gas[:,:2], axis=1)
key = np.where(R < Rcut)[0]
gas_mass = sn_gas.part0.mass
gas_mass[key] = 0.0

# now compute the tidal radius and perform a cut
rt = compute_tidal_radius(R_Sgr, sn, sn_Sgr)
r = np.linalg.norm(sn_Sgr.part1.pos, axis=1)
in_rt = np.where(r < rt)[0]
rs = np.linalg.norm(sn_Sgr.part3.pos, axis=1)
in_rt_star = np.where(rs < rt)[0]

npart[5] = len(in_rt)
masses[5] = sn_Sgr.MassTable[1]

npart[6] = len(in_rt_star)
masses[6] = sn_Sgr.MassTable[3]


#max_id = np.max([np.max(sn_gas.part0.id), np.max(sn.part1.id), np.max(sn.part2.id), np.max(sn.part3.id)])
#TotNTracers, TracerID, ParentID, FluidQuantities = gen_mctracers(gas_mass, sn_gas.part0.id, max_id, Npercell)

#npart[5] = TotNTracers

ics = arepo.ICs('ics.hdf5', npart, masses=masses)

# add MC tracer field
#ics.addField('TracerID', [0, 0, 0, 0, 0, 1], dtype='uint32')
#ics.addField('ParentID', [0, 0, 0, 0, 0, 1], dtype='uint32')
#ics.addField('FluidQuantities', [0, 0, 0, 0, 0, NumFluidQuantities], dtype='<f4')

print(npart)

ics.part0.pos[:] = sn_gas.part0.pos
ics.part0.mass[:] = gas_mass
ics.part0.vel[:] = sn_gas.part0.vel
#ics.part0.id[:]  = sn_gas.part0.id
ics.part0.u[:]   = sn_gas.part0.u

ics.part1.pos[:] = sn.part1.pos + center
ics.part1.vel[:] = sn.part1.vel
#ics.part1.id[:]  = sn.part1.id

ics.part2.pos[:] = sn.part2.pos + center
ics.part2.vel[:] = sn.part2.vel
#ics.part2.id[:]  = sn.part2.id

ics.part3.pos[:] = sn.part3.pos + center
ics.part3.vel[:] = sn.part3.vel
#ics.part3.id[:]  = sn.part3.id

ics.part5.pos[:] = sn_Sgr.part1.pos[in_rt] + center + Sgr_center
ics.part5.vel[:] = sn_Sgr.part1.vel[in_rt] + Sgr_vel

ics.part6.pos[:] = sn_Sgr.part3.pos[in_rt_star] + center + Sgr_center
ics.part6.vel[:] = sn_Sgr.part3.vel[in_rt_star] + Sgr_vel

current_id = 0
for i in range(7):
    if npart[i] > 0:
        part = getattr(ics, 'part'+str(i))
        part.id[:] = np.arange(current_id+1, current_id+1+npart[i])
        current_id = current_id + npart[i]

#max_id = np.max(sn.part3.id)
#ics.part5.id[:] = np.arange(max_id+1, max_id + npart[5] + 1)

ics.write()

