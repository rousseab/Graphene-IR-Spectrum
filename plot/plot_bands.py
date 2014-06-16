
import common
reload(common)
from common import *


#----------------------------------------
# compute dispersion
#----------------------------------------

nkpoints  = 200
kpath     = Path(acell,nkpoints)
ksize     = kpath.list_k.shape[0]

epsilon_k = function_epsilon_k(kpath.list_k)

#----------------------------------------
# plot dispersion
#----------------------------------------

fig = plt.figure(figsize=(14,10))
ax  = fig.add_subplot(111)


ax.plot(kpath.list_x, epsilon_k,'r-',lw=5)
ax.plot(kpath.list_x,-epsilon_k,'b-',lw=5)

#ax.set_title(r'Tight-binding dispersion of Graphene')
ax.set_ylabel(r'$\epsilon_{\bf k}$ (eV)')

ax.set_xticks(kpath.list_xticks)
ax.set_xticklabels(kpath.list_labels)

ax.set_xlim(kpath.list_x[0],kpath.list_x[-1])

ax.grid(True,linestyle='-',color='grey',axis='y',alpha=0.5)
ax.grid(True,linestyle='-',color='black',axis='x',alpha=1.)

fig.subplots_adjust(    left    =       0.15,
                        bottom  =       0.10,
                        right   =       0.95,
                        top     =       0.92,
                        wspace  =       0.50,
                        hspace  =       0.50)

plt.savefig('bands_graphene.pdf')
#plt.show()
