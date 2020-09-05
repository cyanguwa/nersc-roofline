
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

font = { 'size'   : 15}
plt.rc('font', **font)

colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:gray','tab:olive','tab:cyan']
styles = ['o','s','v','^','D',">","<","*","h","H","+","1","2","3","4","8","p","d","|","_",".",","]

markersize = 10
markerwidth = 2
maxchar = 25

def roofline(filename, memRoofs, cmpRoofs, axises, FLOPS, AI, LABELS=None, flag='HBM'):

    if not FLOPS:
        print('FLOPS can not be empty!')
        return
    if max(FLOPS)==0:
        print('FLOPS are all 0s!')
        return

    LABELS = [x[:maxchar] for x in LABELS]

    fig = plt.figure(1,figsize=(10.67,6.6))
    plt.clf()
    ax = fig.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Arithmetic Intensity [FLOPs/Byte]')
    ax.set_ylabel('Performance [GFLOP/sec]')

    nx   = 10000
    xmin = axises[0] #-3 
    xmax = axises[1] #3
    ymin = axises[2] #1
    ymax = axises[3] #5000

    ax.set_xlim(10**xmin, 10**xmax)
    ax.set_ylim(ymin, ymax)

    ixx = int(nx*0.05)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()

    scomp_x_elbow  = []
    scomp_ix_elbow = []
    smem_x_elbow   = []
    smem_ix_elbow  = []

    x = np.logspace(xmin,xmax,nx)
    for roof in cmpRoofs:
        for ix in range(1,nx):
            if float(memRoofs[0][1] * x[ix]) >= roof[1]*1024 and (memRoofs[0][1] * x[ix-1]) < roof[1]*1024:
                scomp_x_elbow.append(x[ix-1])
                scomp_ix_elbow.append(ix-1)
                break

    for roof in memRoofs:
        for ix in range(1,nx):
            if (cmpRoofs[0][1]*1024 <= roof[1] * x[ix] and cmpRoofs[0][1]*1024 > roof[1] * x[ix-1]):
                smem_x_elbow.append(x[ix-1])
                smem_ix_elbow.append(ix-1)
                break

    for i in range(len(cmpRoofs)):
        roof = cmpRoofs[i][1]*1024
        y = np.ones(len(x)) * roof
        ax.plot(x[scomp_ix_elbow[i]:],y[scomp_ix_elbow[i]:],c='k',ls='-',lw='2')

    for i in range(len(memRoofs)):
        roof = memRoofs[i][1]
        y = x * roof
        ax.plot(x[:smem_ix_elbow[i]+1],y[:smem_ix_elbow[i]+1],c='k',ls='-',lw='2')


    marker_handles = []  

    for i in range(len(AI)):
        if AI[i]>0:
#             lab='{:3s} - {:>6.1f} FLOPs/Byte'.format(LABELS[i], AI[i])
            ax.plot(float(AI[i]),float(FLOPS[i]),c=colors[i%10],marker=styles[i],\
                    linestyle='None',ms=markersize,markerfacecolor='none',\
                    markeredgewidth=markerwidth,label=LABELS[i] if LABELS else "unknown")
            marker_handles.append(ax.plot([],[],c=colors[i%10],marker=styles[i],linestyle='None',ms=markersize,\
                    markerfacecolor='none',markeredgewidth=markerwidth,\
                    label='{} - {:.2f} FLOPs/Byte'.format(LABELS[i], AI[i]))[0])

    for roof in cmpRoofs:
        ax.text(x[-ixx],roof[1]*1000,
              roof[0] + ': ' + '{0:.2f}'.format(roof[1]) + ' TFLOP/s',
              horizontalalignment='right',
              verticalalignment='bottom')
    ax.text(x[-ixx],FLOPS[-1],
            '{0:.2f}'.format(FLOPS[-1]/1e3) + ' TFLOP/s',
            horizontalalignment='right',
            verticalalignment='center')

    for roof in memRoofs:
        ang = np.arctan(np.log10(xlim[1]/xlim[0]) / np.log10(ylim[1]/ylim[0])
                                   * fig.get_size_inches()[1]/fig.get_size_inches()[0] )
        if x[ixx]*roof[1] >ymin:
            ax.text(x[ixx],x[ixx]*roof[1]*(1+0.4*np.sin(ang)**2),
              roof[0] + ': ' + '{0:.2f}'.format(float(roof[1])) + ' GB/s',
              horizontalalignment='left',
              verticalalignment='bottom',
              rotation=180/np.pi*ang)
        else:
            ymin_ix_elbow=list()
            ymin_x_elbow=list()
            for ix in range(1,nx):
                if (ymin <= roof[1] * x[ix] and ymin > roof[1] * x[ix-1]):
                    ymin_x_elbow.append(x[ix-1])
                    ymin_ix_elbow.append(ix-1)
                    break
            ax.text(x[ixx+ymin_ix_elbow[0]],x[ixx+ymin_ix_elbow[0]]*roof[1]*(1+0.4*np.sin(ang)**2),
              roof[0] + ': ' + '{0:.2f}'.format(float(roof[1])) + ' GB/s',
              horizontalalignment='left',
              verticalalignment='bottom',
              rotation=180/np.pi*ang)


        
    leg1 = plt.legend(handles = marker_handles,loc='lower right', ncol=1,bbox_to_anchor = (1,0))
    ax.add_artist(leg1)

    plt.savefig('_'.join([filename])+'.png')
#     plt.savefig('_'.join([filename,flag])+'.eps')

    plt.show()



# likwid
memRoofs = [('L1', 1.4*64*2*64), ('L2', 1.4*64*32),  ('HBM', 490), ('DDR', 100)]
cmpRoofs = [('FP64 w/ FMA', 1.4*68*8*2*2/1e3)]
axises=[-3, 4, 1, 5000]
print(memRoofs)
print(cmpRoofs)
LABELS=['L1','L2','HBM','DDR']
flops=5051.923
FLOPS=[flops/10.2243]*4
AI=[flops/6456.799, flops/1387.739, flops/742.8158, flops/0.8883]
print(AI)
print(FLOPS)
roofline('likwid', memRoofs, cmpRoofs, axises, FLOPS, AI, LABELS, 'all')


# sde vtune
memRoofs = [('L1', 1.4*64*2*64), ('L2', 1.4*64*32),  ('HBM', 490), ('DDR', 100)]
cmpRoofs = [('FP64 w/ FMA', 1.4*68*8*2*2/1e3)]
axises=[-3, 4, 1, 5000]
print(memRoofs)
print(cmpRoofs)
LABELS=['L1','L2', 'HBM','DDR']
flops=5839.811
FLOPS=[flops/10.2243]*4
AI=[flops/3795.623, 0, flops/594.562, flops/0.735]
print(AI)
print(FLOPS)
roofline('sde-vtune', memRoofs, cmpRoofs, axises, FLOPS, AI, LABELS, 'all')


# nvprof
memRoofs = [('L1', 54000.), ('L2', 2996.77),  ('HBM', 828.76)]
cmpRoofs = [('FP64', 1312e6*80*32*2/1e12)]
axises=[-3, 4, 100, 15000]
print(memRoofs)
print(cmpRoofs)
LABELS=['L1','L2', 'HBM']
flops=6.3206e+12/1e9
FLOPS=[flops/1.2808043044885498]*3
nvp=[4.8952e+10+328212077, 3.8942e+10+322043992, 1.5013e+10+123354]
byte=np.array(nvp)*32/1e9
AI=flops/byte
print(AI)
print(FLOPS)
roofline('nvprof', memRoofs, cmpRoofs, axises, FLOPS, AI, LABELS, 'all')


# ncu10
memRoofs = [('L1', 54000.), ('L2', 2996.77),  ('HBM', 828.76)]
cmpRoofs = [('FP64', 1.31e9*80*32*2/1e12)]
axises=[-3, 4, 100, 15000]
print(memRoofs)
print(cmpRoofs)
LABELS=['L1','L2', 'HBM']
time=1677685664.15/1.31e9
flops=6.3206e+12/1e9
FLOPS=[flops/time]*3
nvp=[48910700966+322043904,38675595312+322043904+645,15016718515+127241]
byte=np.array(nvp)*32/1e9
AI=flops/byte
print(AI)
print(FLOPS[-1])
roofline('ncu10', memRoofs, cmpRoofs, axises, FLOPS, AI, LABELS, 'all')


# ncu11
memRoofs = [('L1', 54000.), ('L2', 2996.77),  ('HBM', 828.76)]
cmpRoofs = [('FP64', 1.31e9*80*32*2/1e12)]
axises=[-3, 4, 100, 15000]
print(memRoofs)
print(cmpRoofs)
LABELS=['L1','L2', 'HBM']
time=1677853638.88/1.31e9
flops=(164886478848+2363372863488*2+1429016150016)/1e9
FLOPS=[flops/time]*3
nvp=[1.58e12,1.25e12,480.44e9]
byte=np.array(nvp)/1e9
AI=flops/byte
print(AI)
print(FLOPS[-1])
roofline('ncu11', memRoofs, cmpRoofs, axises, FLOPS, AI, LABELS, 'all')


