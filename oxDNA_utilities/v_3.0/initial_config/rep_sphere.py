n_pred=5
n_res=5
l=20
L=50
tot=n_pred*L+n_res*l

r0=22.5 # box is 45
stiff=0.005

trap = open("sphere.txt","w")

for i in range (0,tot):
    trap.write(str(i)+" 0 0 0 "+str(stiff)+" "+str(r0)+" "+str(0)+"\n")

trap.close()
