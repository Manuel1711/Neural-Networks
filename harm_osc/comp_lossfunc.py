import sys
import numpy as np

string = np.array(["00000000000000000000000000000000000000000000000000000000000000000000000000000"]*2)


if len(sys.argv) != 3 :
    print(f"Error - Number of inputs different by 2, but given {len(sys.argv)-1}")
    sys.exit(2)


for i in range(1,len(sys.argv)):
        string[i-1] = sys.argv[i]    

print(string)

rho_nn = np.loadtxt(string[0], unpack='True')
rho_true = np.loadtxt(string[1], unpack='True')

if len(rho_nn[1]) != len(rho_true[1]) :
    print(f"Mismatch between the 2 outputs: {len(rho_nn[1])} vs {len(rho_nn[0])}")
    sys.exit(2)

ftemp = open(file="ttemp.txt", mode='w')
ftemp.write(f"{string[0]}\t\t{string[1]}\n\n")
for i in range(0, len(rho_nn[1])):
    ftemp.write(f"{rho_nn[1][i]}\t\t{rho_true[1][i]}\n\n")
ftemp.close()


print(np.linalg.norm(rho_nn[1] - rho_true[1])**1.)
