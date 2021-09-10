import numpy as np

NStates = 2 # Inclduing ground state
init = 1 # 0, 1, 2
NTraj = 1000

################# PYRAMIDAL ##################

ek = np.zeros((NStates,NTraj))
gamm = 1/3

file01 = open("ek_dist.dat","w")
file02 = open("nk_dist.dat","w")
file01.write("{")
file02.write("{")

for t in range(NTraj):
    while (True):
        rand = np.random.rand()
        ek[init,t] = rand
        rand = np.random.rand()
        if (1 - ek[init,t] > rand):
            break

    for i in range(NStates):
        if (i != init):
            rand = np.random.rand()
            rand = rand * (1 - ek[init,t])
            ek[i,t] = rand
    ek[init] += 1
    
    file01.write("{" + str(ek[0,t]) + ", " + str(ek[1,t]) + "}," + "\n")
    file02.write("{" + str(ek[0,t] - gamm) + ", " + str(ek[1,t] - gamm) + "}," + "\n")

file01.write("}")
file02.write("}")

################# Square ##################

ek = np.zeros((NStates,NTraj))
gamm = (np.sqrt(3) - 1)*0.5

file01 = open("ek_rec_dist.dat","w")
file02 = open("nk_rec_dist.dat","w")
file01.write("{")
file02.write("{")

for t in range(NTraj):
    for i in range(NStates):
        rand = np.random.rand()
        ek[i,t] = 2*gamm*rand
    ek[init] += 1
    
    file01.write("{" + str(ek[0,t]) + ", " + str(ek[1,t]) + "}," + "\n")
    file02.write("{" + str(ek[0,t] - gamm) + ", " + str(ek[1,t] - gamm) + "}," + "\n")

file01.write("}")
file02.write("}")