from sys import argv
from pylab import *

filename = argv[1]
   
data = loadtxt(filename, skiprows=2)
blokSize = data[0,:]
sigma = data[2,:]

plot(blokSize, sigma)
xlabel('Block size')
ylabel('$\sigma$', fontsize=16)
savefig('/home/anette/output/blocking_H2_oppg3.jpg')
show()

