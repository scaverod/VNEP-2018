import numpy as np
import matplotlib.pyplot as plt

#files = [ 'fdensevir1', 'fdensevir4', 'fsparsevir1', 'fsparsevir4', \
#          'fhiervir1', 'fhiervir4', 'ftsvir1', 'ftsvir4']
files = [ 'fdensevir4', 'fsparsevir4', 'fhiervir4', 'ftsvir4']

for filename in files:
    data = np.loadtxt(filename + '.dat', delimiter='&', skiprows=1, dtype=[('s', 'f'), ('v', 'f'), \
            ('btl', 'int'), ('bfi', 'f'), ('bti', 'f'), ('bcost', 'f'), \
            ('ctl', 'int'), ('cfi', 'f'), ('cti', 'f'), ('ccost', 'f')])

    data.sort(order='v')
    print(filename)
    #print(data)
    #print(data['v'], data['s'], data['bti'])

    ratio_s = 2.0 / max(data['s'])
    print(max(data['s']), ratio_s)

    plt.figure(1)
    plt.plot(data['v'] + ratio_s * data['s'], data['bti'], 'b-', label='BP')
    plt.plot(data['v'] + ratio_s * data['s'], data['cti'], 'r--', label='ILP')
    plt.legend(loc='best')
    plt.ylabel('time (s)')
    plt.xlabel('Instance Size')
    plt.xlim(2,20)
    plt.savefig('../graphs/time_' + filename + '.pdf')
    plt.clf()

    plt.figure(2)
    plt.plot(data['v'] + ratio_s * data['s'], data['bcost'], 'b-', label="BP")
    plt.plot(data['v'] + ratio_s * data['s'], data['ccost'], 'r--', label="ILP")
    plt.legend(loc='best')
    plt.ylabel('cost')
    plt.xlabel('Instance Size')
    plt.xlim(2,20)
    plt.savefig('../graphs/cost_' + filename + '.pdf')
    plt.clf()

    plt.figure(3)
    plt.plot(data['v'] + ratio_s * data['s'], 9 - data['btl'], 'b-', label="BP")
    plt.plot(data['v'] + ratio_s * data['s'], 9 - data['ctl'], 'r--', label="ILP")
    plt.ylim(-1,10)
    plt.legend(loc='best')
    plt.ylabel('# Solved Instances')
    plt.xlabel('Instance Size')
    plt.xlim(2,20)
    plt.savefig('../graphs/opt_' + filename + '.pdf')
    plt.clf()

    plt.figure(4)
    plt.plot(data['v'] + ratio_s * data['s'], data['bfi'], 'b-', label="BP")
    plt.plot(data['v'] + ratio_s * data['s'], data['cfi'], 'r--', label="ILP")
    plt.legend(loc='best')
    plt.ylabel('time (s)')
    plt.xlabel('Instance Size')
    plt.xlim(2,20)
    plt.savefig('../graphs/fint_' + filename + '.pdf')
    plt.clf()
