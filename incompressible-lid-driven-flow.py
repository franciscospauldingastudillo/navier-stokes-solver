# This version was last updated: May 28, 2021 @ 6pm
# All run-specific modifications are performed inside NS_solver()

def NS_solver():
    global nx,ny,dx,dy,Lx,Ly,nq,ng,nu,dt,nt,uwall,nqplt,Re,path,restart,trestart
    
    ######## SPATIAL OPTIONS #######
    nx=80
    ny=80
    Lx=1
    Ly=1
    dx=Lx/nx
    dy=Ly/ny
    nq=(nx-1)*ny+nx*(ny-1)
    ng=nx*ny
    nqplt=2*nx*ny
    
    ######## TEMPORAL OPTIONS ######
    dt = 0.0001
    tmin = 0
    tmax = 0.1
    nt = len(np.arange(tmin,tmax+dt,dt))
    uwall=1
    
    ######## VISCOSITY OPTIONS ######
    Re = 400
    nu = 1/Re # implicitly assumes that L=1,uwall=1
    
    ########  RESTART OPTIONS #######
    restart = False
    trestart = 0
    
    ######## DATA/OUTPUT PATH #######
    path = '/home/fspauldinga/MECH-250H/Checkpoint-05/CP5-%dx%d-Re-%d-dt-%1.4f'%(nx,ny,Re,dt)
    if not os.path.exists(path):
        cmd2 = 'mkdir -p %s' % (path)
        os.system(cmd2)
    
    initiate_NS()
    initiate_pp(0)
    #check_matrix()
    integrate_NS()
    print('done.')
    
#############################################################################################################   
def integrate_NS():
    global nx,ny # https://vbsreddy1.medium.com/unboundlocalerror-when-the-variable-has-a-value-in-python-e34e097547d6
    
    uF = np.zeros([nq],dtype='float')   # IC: uF corresponds to fluid at rest
    qnm1 = np.zeros([nq],dtype='float')   
    qn = np.zeros([nq],dtype='float')   # IC: fluid at rest
    gn = np.zeros([ng],dtype='float')   # IC: isostatic pressure field 
    qnp1 = np.zeros([nq],dtype='float')  
    gnp1 = np.zeros([ng],dtype='float')
    
    # Optional: Initialize Simulation with Restart File
    if restart:
        if not os.path.exists(path):
            cmd2 = 'mkdir -p %s' % (path)
            os.system(cmd2)
        t0 = trestart # select restart timestep index to start from
        filename = '%s/CP5-%dx%d-t-%d.pickle'%(path,nx,ny,t0)
        if os.path.exists(filename):
            with open(filename, 'rb') as f:
                qnm1,qn,gn,t0,nx,ny = pickle.load(f)
                print('Initializing %dx%d Re=%d Simulation with Restart File from t='%(nx,ny,Re),t0)
    else:
        t0 = 0
        print('Initializing %dx%d Re=%d Simulation from t='%(nx,ny,Re),t0)
    
    for t in range(t0,nt):
        print('t=%d of nt=%d'%(t,nt))
        # Projection Method: Step 1
        RHS = S(qn[:])+(dt/2)*(3*Ad(qn[:])-Ad(qnm1[:]))+(nu*dt/2)*(BC_L(qn[:])+BC_L(qn[:]))
        #RHS = S(qn[:])+(nu*dt/2)*(BC_L(qn[:])+BC_L(qn[:]))
        guess = qn[:] # current velocity
        uF[:] = conjgrad(RHS,guess,imax=100,eps=10**(-6),opt='mom-eq')
        # Projection Method: Step 2 
        RHS = (1/dt)*div(uF[:])+(1/dt)*BC_D()
        guess = gn[:]
        gnp1[:] = conjgrad(-RHS,guess,imax=10**4,eps=10**(-6),opt='pp-eq') # negative definite to positive definite
        gnp1[0]=pp[0,0] # pin pressure
        
        # Projection Method: Step 3
        RHS = uF[:]-dt*Rinv(grad(gnp1[:]))
        qnp1[:] = RHS
        
        # Prepare fields for next time step
        qnm1 = qn
        qn = qnp1
        gn = gnp1
        qnp1 = np.zeros([nq],dtype='float')  
        gnp1 = np.zeros([ng],dtype='float')
        uF   = np.zeros([nq],dtype='float')
        
        if errimax==1:
            # CGM solver failed on Pressure Poisson
            # Abort the Simulation
            print('CGM solver failed on Pressure Poisson: Aborting Run at time=',t*dt)
            break

        # INSERT PLOTTING ROUTINE HERE
        if ((t%50==0 or t==nt-1) and nx==ny and t!=t0):
            ## The conditional statement (t!=t0) prevents an overwrite of the original pickle file/plot.
            # Save Pickle File at Regular Intervals
            t0 = t
            filename = '%s/CP5-%dx%d-t-%d.pickle'%(path,nx,ny,t0)
            with open(filename, 'wb') as f:
                pickle.dump([qnm1,qn,gn,t0,nx,ny], f)
                print('Saving Restart File at t=',t0)
            
            qplt = np.zeros([nqplt],dtype='float')
            for i in range(0,nx-1):
                for j in range(0,ny):
                    qplt[uplt[i,j]] = qn[u[i,j]]
            for i in range(0,nx):
                for j in range(0,ny-1):
                    qplt[vplt[i,j]] = qn[v[i,j]]
            i=nx-1
            for j in range(0,ny):
                qplt[uplt[i,j]] = 0
            j=ny-1
            for i in range(0,nx):
                qplt[vplt[i,j]] = 0

            xg = np.empty([nx],dtype='float')
            yg = np.empty([ny],dtype='float')
            for i in range(0,nx):
                xg[i] = i*dx+dx/2
            for j in range(0,ny):
                yg[j] = j*dy+dy/2

            plt.figure(7)
            plt.quiver(xg,yg,qplt[uplt.T],qplt[vplt.T],color='r',linewidth=0.5,headwidth=3)
            ax = plt.gca()               
            ax.set_ylabel('ny')
            ax.set_xlabel('nx') 
            plt.show()
            figname = '%s/CP5-uv-%dx%d-T-%f.png'%(path,nx,ny,t*dt)
            plt.savefig(figname)

            xg = np.empty([nx],dtype='float')
            yg = np.empty([ny],dtype='float')
            M = np.empty([nx,ny],dtype='float')
            for i in range(0,nx):
                xg[i] = i*dx+dx/2
            for j in range(0,ny):
                yg[j] = j*dy+dy/2
            for i in range(0,nx): 
                for j in range(0,ny):
                    if (i==0 and j==0):
                        M[i,j] = pp[0,0]
                    else:
                        M[i,j] = gn[p[j,i]]

            plt.figure(8)
            im = plt.contourf(xg,yg,M,color='r',linewidth=0.5,headwidth=3)
            ax = plt.gca()               
            ax.set_ylabel('ny')
            ax.set_xlabel('nx')
            plt.colorbar(im)
            plt.show()
            figname = '%s/CP5-p-%dx%d-T-%f.png'%(path,nx,ny,t*dt)
            plt.savefig(figname)

            # Plot U and V Separately
            plt.figure(12)
            im = plt.imshow(qplt[uplt.T])
            ax = plt.gca()               
            ax.set_ylabel('ny')
            ax.set_xlabel('nx')
            ax.invert_yaxis()
            plt.colorbar(im)
            plt.show()
            figname = '%s/CP5-u-%dx%d-T-%f.png'%(path,nx,ny,t*dt)
            plt.savefig(figname)

            plt.figure(13)
            im = plt.imshow(qplt[vplt.T])
            ax = plt.gca()               
            ax.set_ylabel('ny')
            ax.set_xlabel('nx')
            ax.invert_yaxis()
            plt.colorbar(im)
            plt.show()
            figname = '%s/CP5-v-%dx%d-T-%f.png'%(path,nx,ny,t*dt)
            plt.savefig(figname)
        
            plt.close('all')
        
        # Optional: put diagnostic information in an output file
        if 1: 
            # https://cmdlinetips.com/2012/09/three-ways-to-write-text-to-a-file-in-python/
            # https://stackoverflow.com/questions/415511/how-to-get-the-current-time-in-python
            outfile = open('%s/OUTPUT.txt'%(path), 'w')
            str1 = 't=%d of nt=%d'%(t,nt)
            str2 = str(datetime.now())
            lineList = [str1,str2]
            for line in lineList:
                # write each line to output file
                outfile.write(line)
                outfile.write("\n")
            outfile.close()
      
     
            
#############################################################################################################   
def R(q):
    # Input(q): size velocity
    # Output(r): size velocity
    r = np.zeros([nq],dtype='float')
    r = q - (dt*nu/2)*lap(q)
    return r

#############################################################################################################   
def S(q):
    # Input(q): size velocity
    # Output(s): size velocity
    s = np.zeros([nq],dtype='float')
    s = q + (dt*nu/2)*lap(q)
    return s

#############################################################################################################   
def Rinv(q):
    # Input(q): size velocity
    # Output(rinv): size velocity
    rinv = np.zeros([nq],dtype='float')
    c = (dt*nu/2)
    Lq = lap(q)
    rinv = q + c*Lq + c**2*lap(Lq)
    return rinv

#############################################################################################################   
def Ax_new(xm1):
    # Input(xm1): size pressure minus one
    # Output(Axx): size pressure minus one (not solving for pinned)
    # Here, we evaluate A for the conjugate gradient solver of the pp-equation.
    # -LHS = -RHS
    # This ensures that A is positive definite.
    Bxx = np.zeros([ng],dtype='float')
    x = np.zeros([ng],dtype='float')
    x[1:-1] = xm1
    Bxx = -div(Rinv(grad(x)))
    Axx = np.zeros([ng-1],dtype='float')
    Axx = Bxx[1:-1]
    return Axx

#############################################################################################################   
def Ax(x):
    # Input(xm1): size pressure minus one
    # Output(Axx): size pressure minus one (not solving for pinned)
    # Here, we evaluate A for the conjugate gradient solver of the pp-equation.
    # -LHS = -RHS
    # This ensures that A is positive definite.
    Axx = -div(Rinv(grad(x)))
    return Axx

#############################################################################################################   
def conjgrad_new(b,x,imax,eps,opt):
    # imax: maximum number of iterations
    # eps: error tolerance
    # b: RHS of mom-/pp- equation
    # x: initial guess of solution to mom-/pp- equation
    global errimax
    errimax=0
    delta_new = np.dtype('float64')
    delta_old = np.dtype('float64')
    delta_0 = np.dtype('float64')
    beta = np.dtype('float64')
    alpha = np.dtype('float64')
    
    if opt=='pp-eq':
        # Inputs(x&b) - size pressure
        # NOTE: Here, the initial guess for x must be p_{n+1}
        # Criteria for end of loop is imax and r
        r = np.zeros([ng-1],dtype='float64') # residual
        d = np.zeros([ng-1],dtype='float64') # search direction
        q = np.zeros([ng-1],dtype='float64') # 'A' operator acting on search direction
        i = 0
        r = b[1:-1]-Ax(x[1:-1])            # r = b-Ax = f'(x). Here, A is the operator on pressure at t=n+1
        d = r                              # initial search direction is 'downhill'
        delta_new = np.matmul(r,r)         # initial residual magnitude
        delta_0 = delta_new 
        while (i<imax and delta_new>(delta_0*eps**2)):
            q = Ax(d)                           # q = Ad
            alpha = delta_new/(np.matmul(d,q))  # compute correction to x
            x[1:-1] = x[1:-1] + alpha*d         # update x guess
            if (i%50==0):
                r = b[1:-1]-Ax(x[1:-1]) 
            else:
                r = r-alpha*q                   # compute new residual
            delta_old = delta_new
            delta_new = np.matmul(r,r)          # compute magnitude of new residual
            beta = delta_new/delta_old          # compute correction to search direction
            d = r + beta*d                      # compute new search direction
            i=i+1
    elif opt=='mom-eq':
        # Inputs(x,b) - size velocity
        # NOTE: Here, the initial guess for x must be u_F
        r = np.zeros([nq],dtype='float64') # residual
        d = np.zeros([nq],dtype='float64') # search direction
        q = np.zeros([nq],dtype='float64') # 'A' operator acting on search direction
        i = 0
        r = b-R(x)                       # r = b-Ax = f'(x). Here, A is the operator on u_F at t=n+F<n+1
        d = r                            # initial search direction is 'downhill'
        delta_new = np.matmul(r,r)       # initial residual magnitude
        delta_0 = delta_new 
        while (i<imax and delta_new>(delta_0*eps**2)):
            q = R(d)                            # q = Ad
            alpha = delta_new/(np.matmul(d,q))  # compute correction to x
            x = x + alpha*d                     # update x guess
            if (i%50==0):
                r = b-R(x)
            else:
                r = r - alpha*q                 # compute new residual
            delta_old = delta_new
            delta_new = np.matmul(r,r)          # compute magnitude of new residual
            beta = delta_new/delta_old          # compute correction to search direction
            d = r + beta*d                      # compute new search direction
            i=i+1
    else:
        print('conjgrad: incorrect option specified')
        time.sleep(30)
    
    if (i<imax):
        print('conjgrad success: %s, criteria met at i=%d'%(opt,i-1))
    elif (i==imax):
            print('conjgrad failure: %s, i=imax=%d, change initial guess (x) to %s'%(opt,imax,opt))
            errimax=1
    return x

#############################################################################################################   
def conjgrad(b,x,imax,eps,opt):
    ##### NOTE: This is the old verison of the pressure poisson solver!!! It doesn't remove p[0] from solver.
    
    # imax: maximum number of iterations
    # eps: error tolerance
    # b: RHS of mom-/pp- equation
    # x: initial guess of solution to mom-/pp- equation
    
    global errimax
    errimax=0
    delta_new = np.dtype('float64')
    delta_old = np.dtype('float64')
    delta_0 = np.dtype('float64')
    beta = np.dtype('float64')
    alpha = np.dtype('float64')
    
    if opt=='pp-eq':
        # Inputs(x&b) - size pressure
        # NOTE: Here, the initial guess for x must be p_{n+1}
        r = np.zeros([ng],dtype='float64') # residual
        d = np.zeros([ng],dtype='float64') # search direction
        q = np.zeros([ng],dtype='float64') # 'A' operator acting on search direction
        i = 0
        r = b-Ax(x)                      # r = b-Ax = f'(x). Here, A is the operator on pressure at t=n+1
        d = r                            # initial search direction is 'downhill'
        delta_new = np.matmul(r,r)       # initial residual magnitude
        delta_0 = delta_new 
        while (i<imax and delta_new>(delta_0*eps**2)):
            q = Ax(d)                           # q = Ad
            alpha = delta_new/(np.matmul(d,q))  # compute correction to x
            x = x + alpha*d                     # update x guess
            if (i%50==0):
                r = b-Ax(x) 
            else:
                r = r-alpha*q                   # compute new residual
            delta_old = delta_new
            delta_new = np.matmul(r,r)          # compute magnitude of new residual
            beta = delta_new/delta_old          # compute correction to search direction
            d = r + beta*d                      # compute new search direction
            i=i+1
    elif opt=='mom-eq':
        # Inputs(x,b) - size velocity
        # NOTE: Here, the initial guess for x must be u_F
        r = np.zeros([nq],dtype='float64') # residual
        d = np.zeros([nq],dtype='float64') # search direction
        q = np.zeros([nq],dtype='float64') # 'A' operator acting on search direction
        i = 0
        r = b-R(x)                       # r = b-Ax = f'(x). Here, A is the operator on u_F at t=n+F<n+1
        d = r                            # initial search direction is 'downhill'
        delta_new = np.matmul(r,r)       # initial residual magnitude
        delta_0 = delta_new 
        while (i<imax and delta_new>(delta_0*eps**2)):
            q = R(d)                            # q = Ad
            alpha = delta_new/(np.matmul(d,q))  # compute correction to x
            x = x + alpha*d                     # update x guess
            if (i%50==0):
                r = b-R(x)
            else:
                r = r - alpha*q                 # compute new residual
            delta_old = delta_new
            delta_new = np.matmul(r,r)          # compute magnitude of new residual
            beta = delta_new/delta_old          # compute correction to search direction
            d = r + beta*d                      # compute new search direction
            i=i+1
    else:
        print('conjgrad: incorrect option specified')
        time.sleep(30)
    
    if (i<imax):
        print('conjgrad success: %s, criteria met at i=%d'%(opt,i-1))
    elif (i==imax):
            print('conjgrad failure: %s, i=imax=%d, change initial guess (x) to %s'%(opt,imax,opt))
            errimax=1
    return x
    
#############################################################################################################   
def initiate_NS():
    global u,v,p,w,uplt,vplt
    
    u = np.empty([nx-1,ny],dtype='int')
    v = np.empty([nx,ny-1],dtype='int')
    p = np.empty([nx,ny],dtype='int')
    w = np.empty([nx,ny],dtype='int') # vorticity @ cell corners
    
    #plt.close('all')
    #plt.figure(1)
    #ax = plt.gca()
    
    # create pointers u,v
    index=-1
    for i in range(0,nx-1):
        for j in range(0,ny):
            index=index+1
            u[i,j]=index
            #ax.scatter(i*dx+dx,j*dy+dy/2,marker='>',facecolors='none',edgecolors='b')
    for i in range(0,nx):
        for j in range(0,ny-1):
            index = index+1
            v[i,j]=index
            #ax.scatter(i*dx+dx/2,j*dy+dy,marker='^',facecolors='none',edgecolors='r')
    if index!=((nx-1)*ny+nx*(ny-1)-1):
        print('error: wrong velocity size')
    # create pointers for p
    index=-1
    for i in range(0,nx):
        for j in range(0,ny):
            index=index+1
            if (i==0 and j==0):
                #skip pinned pressure
                p[i,j]=10**(9)
            else:
                p[i,j]=index
                #ax.scatter(i*dx+dx/2,j*dy+dy/2,marker='o',facecolors='none',edgecolors='g') 
    if index!=(nx*ny-1):
        print('error: wrong pressure size')
    # create pointers for w (vorticity)
    index=-1
    for i in range(0,nx):
        for j in range(0,ny):
            index=index+1
            w[i,j]=index
            #ax.scatter(i*dx+dx,j*dy+dy,marker='s',color='k')
    if index!=(nx*ny-1):
        print('error: wrong vorticity size')    
        
    # Visualize Grid #
    #ax.set_ylabel('y-axis')
    #ax.set_xlabel('x-axis') 
    #xticks = np.arange(0,Lx,dx)
    #yticks = np.arange(0,Ly,dy)
    #ax.set_xticks(xticks)
    #ax.set_yticks(yticks)
    #ax.set_xlim(0,Lx)
    #ax.set_ylim(0,Ly)
    #ax.grid(which='both')
    #plt.show()
    
    uplt = np.empty([nx,ny],dtype='int')
    vplt = np.empty([nx,ny],dtype='int')
    # create pointers uplt,vplt
    index=-1
    for i in range(0,nx):
        for j in range(0,ny):
            index=index+1
            uplt[i,j]=index
    for i in range(0,nx):
        for j in range(0,ny):
            index = index+1
            vplt[i,j]=index
    if index!=(2*nx*ny-1):
        print('error: wrong uplt and vplt size')
    

#############################################################################################################
    
def curlz(q):
    ## This routine solves for the z-component of curl at top-right corner of a cell.
    ## It is called as part of the plotting routine. It is not verified for use in time-integration. 
    # The top and right walls introduce BC and GC. 
    # The vorticity @ the top-right wall is determined entirely from BC and GC.
    
    # Input(q) ~ size velocity 
    # Output(ddq) ~ size velocity
    cq = np.zeros([nq],dtype='float')
    
    # Case 1: bottom-left corner
    i = 0
    j = 0
    cq[w[i,j]] = (
                 (-q[v[i,j]]+q[v[i+1,j]])/dx -
                 (-q[u[i,j]]+q[u[i,j+1]])/dy 
                 )
    
    # Case 2: bottom-right corner
    i = nx-1
    j = 0
    cq[w[i,j]] = (
                 (-q[v[i,j]]  +gcv(False,q[v[i,j]]))/dx - # INSERT GC: +q[v[i+1,j]]/dx
                 (-bcu(i,j)+bcu(i,j+1))/dy             # INSERT BC: -q[u[i,j]]/dy,+q[u[i,j+1]]/dy
                 )
    
    # Case 3: top-left corner
    i = 0
    j = ny-1
    cq[w[i,j]] = (
                 (-bcv(i,j)+bcv(i+1,j))/dx -            # INSERT BC: -q[v[i,j]]/dx,+q[v[i+1,j]]/dx
                 (-q[u[i,j]]  +gcu(True,q[u[i,j]]))/dy     # INSERT GC: +q[u[i,j+1]]/dy
                 )
   
    # Case 4: top-right corner ***
    i = nx-1
    j = ny-1
    cq[w[i,j]] = (
                 (-bcv(i,j)+gcv(True,bcv(i,j)))/dx -   # INSERT BC: -q[v[i,j]]/dx=0, GC: +q[v[i+1,j]]/dx=0
                 (-bcu(i,j)+gcu(True,bcu(i,j)))/dy     # INSERT BC: -q[u[i,j]]/dy=0, GC: +q[u[i,j+1]]/dy=2uwall
                 )

    # Case 5: bottom wall
    j = 0
    for i in range(1,nx-1):
        cq[w[i,j]] = (
                     (-q[v[i,j]]+q[v[i+1,j]])/dx -
                     (-q[u[i,j]]+q[u[i,j+1]])/dy 
                     )

    # Case 6: top wall
    j = ny-1
    for i in range(1,nx-1):
        cq[w[i,j]] = (
                     (-bcv(i,j)+bcv(i+1,j))/dx -          # INSERT BC: -q[v[i,j]]/dx,+q[v[i+1,j]]/dx
                     (-q[u[i,j]]  +gcu(True,q[u[i,j]]))/dy   # INSERT GC: +q[u[i,j+1]]/dy
                     )
    
    # Case 7: left wall
    i = 0
    for j in range(1,ny-1):
        cq[w[i,j]] = (
                     (-q[v[i,j]]+q[v[i+1,j]])/dx -
                     (-q[u[i,j]]+q[u[i,j+1]])/dy 
                     )

    # Case 8: right wall
    i = nx-1
    for j in range(1,ny-1):
        cq[w[i,j]] = (
                     (-q[v[i,j]]  +gcv(False,q[v[i,j]]))/dx - # INSERT GC: +q[v[i+1,j]]/dx
                     (-bcu(i,j)+bcu(i,j+1))/dy             # INSERT BC: -q[u[i,j]]/dy,+q[u[i,j+1]]/dy
                     )
        
    # Case 9: interior
    for i in range(1,nx-1):
        for j in range(1,ny-1):
            cq[w[i,j]] = (
                         (-q[v[i,j]]+q[v[i+1,j]])/dx -
                         (-q[u[i,j]]+q[u[i,j+1]])/dy 
                         )
    return cq   
    
    
#############################################################################################################
    
def div(q):
    # Input(q) ~ size velocity 
    # Output(g) ~ size pressure
    g = np.zeros([ng],dtype='float')
    # Note: this routine solves for 'g' at cell centers
    
    # Case 1: bottom-left corner
    i = 0
    j = 0
    # p(0,0) is pinned and not a part of div()
    # -- do not evaluate.
    
    # Case 2: bottom-right corner
    i = nx-1
    j = 0
    g[p[i,j]] = (
                (-q[u[i-1,j]]          )/dx +     # INSERT BC: +q[u[i,j]]/dx
                (            +q[v[i,j]])/dy       # INSERT BC: -q[v[i,j-1]]/dy
                )
    
    # Case 3: top-left corner
    i = 0
    j = ny-1
    g[p[i,j]] = (
                (            +q[u[i,j]])/dx +     # INSERT BC: -q[u[i-1,j]]/dx
                (-q[v[i,j-1]]          )/dy       # INSERT BC: +q[v[i,j]]/dy
                )
    
    # Case 4: top-right corner
    i = nx-1
    j = ny-1
    g[p[i,j]] = (
                (-q[u[i-1,j]]          )/dx +     # INSERT BC: +q[u[i,j]]/dx
                (-q[v[i,j-1]]          )/dy       # INSERT BC: +q[v[i,j]]/dy
                )
    
    # Case 5: bottom wall
    j = 0
    i = np.arange(1,nx-1)
    #for i in range(1,nx-1):
    g[p[i,j]] = (
                (-q[u[i-1,j]]+q[u[i,j]])/dx +
                (            +q[v[i,j]])/dy   # INSERT BC: -q[v[i,j-1]]/dy
                )

    # Case 6: top wall
    j = ny-1
    i = np.arange(1,nx-1)
    #for i in range(1,nx-1):
    g[p[i,j]] = (
                (-q[u[i-1,j]]+q[u[i,j]])/dx +
                (-q[v[i,j-1]]          )/dy   # INSERT BC: +q[v[i,j]]/dy
                )
    
    # Case 7: left wall
    i = 0
    j = np.arange(1,ny-1)
    #for j in range(1,ny-1):
    g[p[i,j]] = (
                (            +q[u[i,j]])/dx + # INSERT BC: -q[u[i-1,j]]/dx
                (-q[v[i,j-1]]+q[v[i,j]])/dy 
                )
    
    # Case 8: right wall
    i = nx-1
    j = np.arange(1,ny-1)
    #for j in range(1,ny-1):
    g[p[i,j]] = (
                (-q[u[i-1,j]]          )/dx + # INSERT BC: +q[u[i,j]]/dx
                (-q[v[i,j-1]]+q[v[i,j]])/dy 
                )
    
    # Case 9: interior
    i = np.arange(1,nx-1)
    j = np.arange(1,ny-1)
    #for i in range(1,nx-1):
        #for j in range(1,ny-1):
    g[p[1:nx-1,1:ny-1]] = np.add(
                                (-q[u[i-1,1:ny-1]]+q[u[i,1:ny-1]])/dx,
                                (-q[v[1:nx-1,j-1]]+q[v[1:nx-1,j]])/dy 
                                )
        

    
    return g

#############################################################################################################
def BC_D():
    # Output(g) ~ size pressure
    g = np.zeros([ng],dtype='float')
    
    ###### From Checkpoint 02 ######
    # Define a vector u = [sin(x),sin(y)]
    # Then, div(u) = cos(x) + cos(y)
    
    ## GENERAL FORM: locations of velocity edgepoints
    # x-dir
    #xu = i*dx+dx
    #u = sin(xu)
    # y-dir
    #yv = j*dy+dy 
    #v = sin(yv)
    
    ####
    # NOTE: This is the exception because I am calculating the exact derivative right here
    # Case 1: bottom-left corner
    #i = 0
    #j = 0
    #xu = i*dx+dx
    #yv = j*dy+dy 
    #g[0] = np.cos(xu)+np.cos(yv) # this is the pinned value
    ####
    ###### ^ From Checkpoint 02 ^ ######
    
    # Case 1: bottom-left corner
    i = 0
    j = 0
    # p(0,0) is pinned and not a part of BC_D()
    # -- do not evaluate.
    
    # Case 2: bottom-right corner
    i = nx-1
    j = 0
    g[p[i,j]] = (
                ( bcu(i,j)/dx) +   # INSERT BC: +q[u[i,j]]/dx
                (-bcv(i,j-1)/dy)   # INSERT BC: -q[v[i,j-1]]/dy
                )
    
    # Case 3: top-left corner
    i = 0
    j = ny-1
    g[p[i,j]] = (
                (-bcu(i-1,j)/dx) +   # INSERT BC: -q[u[i-1,j]]/dx
                ( bcv(i,j)/dy)       # INSERT BC: +q[v[i,j]]/dy
                )
    
    # Case 4: top-right corner
    i = nx-1
    j = ny-1
    g[p[i,j]] = (
                (bcu(i,j)/dx) +     # INSERT BC: +q[u[i,j]]/dx
                (bcv(i,j)/dy)       # INSERT BC: +q[v[i,j]]/dy
                )
    
    # Case 5: bottom wall
    j = 0
    i = np.arange(1,nx-1)
    #for i in range(1,nx-1):
    g[p[i,j]] = (
                (-bcv(i,j-1)/dy)   # INSERT BC: -q[v[i,j-1]]/dy
                )

    # Case 6: top wall
    j = ny-1
    i = np.arange(1,nx-1)
    #for i in range(1,nx-1):
    g[p[i,j]] = (
                (bcv(i,j)/dy)   # INSERT BC: +q[v[i,j]]/dy
                )
    
    # Case 7: left wall
    i = 0
    j = np.arange(1,ny-1)
    #for j in range(1,ny-1):
    g[p[i,j]] = (
                (-bcu(i-1,j)/dx)  # INSERT BC: -q[u[i-1,j]]/dx
                )
    
    # Case 8: right wall
    i = nx-1
    j = np.arange(1,ny-1)
    #for j in range(1,ny-1):
    g[p[i,j]] = (
                (bcu(i,j)/dx)  # INSERT BC: +q[u[i,j]]/dx
                )
    
    # Case 9: interior (no BC to evalutate)
    
    return g   

#############################################################################################################
def initiate_pp(val):
    # This subroutine intitiates the pinned pressure as a GLOBAL 2D array.
    # Note: pp is used by grad()
    # 05/15/21: pp[bottom-left corner] = 1
    global pp
    pp = np.zeros([nx,ny],dtype='float')
    pp[0,0] = val

#############################################################################################################
#def pp(i,j):
    # NOTE: This subroutine MUST be consistent with the function used in grad()
    # 05/08/21: f = sin(x) + sin(y)
    # NOTE: This routine was used in checkpoint-02. It will be phased out.
    #xg = i*dx+dx/2
    #yg = j*dy+dy/2
    #pin = np.sin(xg)+np.sin(yg)    
    #return pin
    
#############################################################################################################

def grad(g):
    # Input(g) ~ size pressure
    # Output(q) ~ size velocity 
    # NOTE: pp(i,j) should either be a global array or a function (we don't want multiple calls)
    q = np.zeros([nq],dtype='float')
        
    # Case 1: interior (x-dir)
    
    i=0
    j=0
    # pinned pressure
    q[u[i,j]] = (
                (-pp[i,j]+g[p[i+1,j]])/dx # INSERT BC: -g[p[i,j]]
                )
    
    i=0
    j=np.arange(1,ny)
    #for j in range(1,ny):
    q[u[i,j]] = (
                (-g[p[i,j]]+g[p[i+1,j]])/dx 
                )
    j=0
    i=np.arange(1,nx-1)
    #for i in range(1,nx-1):
    q[u[i,j]] = (
                (-g[p[i,j]]+g[p[i+1,j]])/dx 
                )
    i=np.arange(1,nx-1)
    j=np.arange(1,ny)
    #for i in range(1,nx-1):
        #for j in range(1,ny):
    q[u[1:nx-1,1:ny]] = (
                (-g[p[i,1:ny]]+g[p[i+1,1:ny]])/dx 
                )
    
    # Case 2: interior (y-dir)
    i=0
    j=0
    # pinned pressure
    q[v[i,j]] = (
                (-pp[i,j]+g[p[i,j+1]])/dy  # INSERT BC: -g[p[i,j]]
                )
    
    i=0
    j=np.arange(1,ny-1)
    #for j in range(1,ny-1):
    q[v[i,j]] = (
                (-g[p[i,j]]+g[p[i,j+1]])/dy 
                )
    j=0
    i=np.arange(1,nx)
    #for i in range(1,nx):
    q[v[i,j]] = (
                (-g[p[i,j]]+g[p[i,j+1]])/dy 
                )
    i=np.arange(1,nx)
    j=np.arange(1,ny-1)
    #for i in range(1,nx):
        #for j in range(1,ny-1):
    q[v[1:nx,1:ny-1]] = (
                (-g[p[1:nx,j]]+g[p[1:nx,j+1]])/dy 
                )
    return q

#############################################################################################################

def lap(q):
    # Input(q) ~ size velocity
    # Output(ddq) ~ size velocity
    ddq = np.zeros([nq],dtype='float')
    
    ## X-DIR
    # Case 1: bottom-left corner
    i=0
    j=0
    ddq[u[i,j]] = (
                  (           -2*q[u[i,j]]+q[u[i+1,j]])/dx**2 + # INSERT BC: q[u[i-1,j]]/dx**2
                  (           -2*q[u[i,j]]+q[u[i,j+1]])/dy**2   # INSERT BC: q[u[i,j-1]]/dx**2
                  )
    # Case 2: bottom-right corner
    i=nx-2
    j=0
    ddq[u[i,j]] = (
                  (q[u[i-1,j]]-2*q[u[i,j]]            )/dx**2 + # INSERT BC: q[u[i+1,j]]/dx**2
                  (           -2*q[u[i,j]]+q[u[i,j+1]])/dy**2   # INSERT BC: q[u[i,j-1]]/dy**2
                  )
    
    # Case 5: bottom wall (excluding corners)
    j=0
    i=np.arange(1,nx-2)
    #for i in range(1,nx-2):
    ddq[u[i,j]] = (
                  (q[u[i-1,j]]-2*q[u[i,j]]+q[u[i+1,j]])/dx**2 +
                  (           -2*q[u[i,j]]+q[u[i,j+1]])/dy**2   # INSERT BC: q[u[i,j-1]]/dy**2
                  )     
    
    # Case 3: top-left corner
    i=0
    j=ny-1
    ddq[u[i,j]] = (
                  (           -2*q[u[i,j]]+q[u[i+1,j]])/dx**2 + # INSERT BC: q[u[i-1,j]]/dx**2
                  (q[u[i,j-1]]-2*q[u[i,j]]            )/dy**2   # INSERT BC: q[u[i,j+1]]/dy**2
                  )
    # Case 4: top-right corner
    i=nx-2
    j=ny-1
    ddq[u[i,j]] = (
                  (q[u[i-1,j]]-2*q[u[i,j]]            )/dx**2 + # INSERT BC: q[u[i+1,j]]/dx**2
                  (q[u[i,j-1]]-2*q[u[i,j]]            )/dy**2   # INSERT BC: q[u[i,j+1]]/dy**2
                  )
    # Case 6: top wall (excluding corners)
    j=ny-1
    i=np.arange(1,nx-2)
    #for i in range(1,nx-2):
    ddq[u[i,j]] = (
                  (q[u[i-1,j]]-2*q[u[i,j]]+q[u[i+1,j]])/dx**2 +
                  (q[u[i,j-1]]-2*q[u[i,j]]            )/dy**2   # INSERT BC: q[u[i,j+1]]/dy**2
                  )
    # Case 7: left wall (excluding corners)
    i=0
    j=np.arange(1,ny-1)
    #for j in range(1,ny-1):
    ddq[u[i,j]] = (
                  (           -2*q[u[i,j]]+q[u[i+1,j]])/dx**2 +  # INSERT BC: q[u[i-1,j]]/dx**2
                  (q[u[i,j-1]]-2*q[u[i,j]]+q[u[i,j+1]])/dy**2 
                  )  
    
    # Case 8: right wall (excluding corners)
    i=nx-2
    j=np.arange(1,ny-1)
    #for j in range(1,ny-1):
    ddq[u[i,j]] = (
                  (q[u[i-1,j]]-2*q[u[i,j]]            )/dx**2 +  # INSERT BC: q[u[i+1,j]]/dx**2
                  (q[u[i,j-1]]-2*q[u[i,j]]+q[u[i,j+1]])/dy**2 
                  )
    # Case 9: interior
    # https://stackoverflow.com/questions/53551421/using-two-indexing-arrays-on-a-numpy-array
    i=np.arange(1,nx-2)
    j=np.arange(1,ny-1)
    #for i in range(1,nx-2):
        #for j in range(1,ny-1):
    ddq[u[1:nx-2,1:ny-1]] = np.add(
                                  (q[u[i-1,1:ny-1]]-2*q[u[i,1:ny-1]]+q[u[i+1,1:ny-1]])/dx**2,
                                  (q[u[1:nx-2,j-1]]-2*q[u[1:nx-2,j]]+q[u[1:nx-2,j+1]])/dy**2 
                                  )

    ## Y-DIR       
    # Case 1: bottom-left corner
    i=0
    j=0
    ddq[v[i,j]] = (
                  (           -2*q[v[i,j]]+q[v[i+1,j]])/dx**2 + # INSERT BC: q[v[i-1,j]]/dx**2
                  (           -2*q[v[i,j]]+q[v[i,j+1]])/dy**2   # INSERT BC: q[v[i,j-1]]/dx**2
                  )
    # Case 2: bottom-right corner
    i=nx-1
    j=0
    ddq[v[i,j]] = (
                  (q[v[i-1,j]]-2*q[v[i,j]]            )/dx**2 + # INSERT BC: q[v[i+1,j]]/dx**2
                  (           -2*q[v[i,j]]+q[v[i,j+1]])/dy**2   # INSERT BC: q[v[i,j-1]]/dy**2
                  )
    # Case 5: bottom wall (excluding corners)
    j=0
    i=np.arange(1,nx-1)
    #for i in range(1,nx-1):
    ddq[v[i,j]] = (
                  (q[v[i-1,j]]-2*q[v[i,j]]+q[v[i+1,j]])/dx**2 +
                  (           -2*q[v[i,j]]+q[v[i,j+1]])/dy**2   # INSERT BC: q[v[i,j-1]]/dy**2
                  )
        
    # Case 3: top-left corner
    i=0
    j=ny-2
    ddq[v[i,j]] = (
                  (           -2*q[v[i,j]]+q[v[i+1,j]])/dx**2 + # INSERT BC: q[v[i-1,j]]/dx**2
                  (q[v[i,j-1]]-2*q[v[i,j]]            )/dy**2   # INSERT BC: q[v[i,j+1]]/dy**2
                  )
    # Case 4: top-right corner
    i=nx-1
    j=ny-2
    ddq[v[i,j]] = (
                  (q[v[i-1,j]]-2*q[v[i,j]]            )/dx**2 + # INSERT BC: q[v[i+1,j]]/dx**2
                  (q[v[i,j-1]]-2*q[v[i,j]]            )/dy**2   # INSERT BC: q[v[i,j+1]]/dy**2
                  )
    # Case 6: top wall (excluding corners)
    j=ny-2
    i=np.arange(1,nx-1)
    #for i in range(1,nx-1):
    ddq[v[i,j]] = (
                  (q[v[i-1,j]]-2*q[v[i,j]]+q[v[i+1,j]])/dx**2 +
                  (q[v[i,j-1]]-2*q[v[i,j]]            )/dy**2   # INSERT BC: q[v[i,j+1]]/dy**2
                  ) 
    
    # Case 7: left wall (excluding corners)
    i=0
    j=np.arange(1,ny-2)
    #for j in range(1,ny-2):
    ddq[v[i,j]] = (
                  (           -2*q[v[i,j]]+q[v[i+1,j]])/dx**2 +  # INSERT BC: q[v[i-1,j]]/dx**2
                  (q[v[i,j-1]]-2*q[v[i,j]]+q[v[i,j+1]])/dy**2 
                  )
    # Case 8: right wall (excluding corners)
    i=nx-1
    j=np.arange(1,ny-2)
    #for j in range(1,ny-2):
    ddq[v[i,j]] = (
                  (q[v[i-1,j]]-2*q[v[i,j]]            )/dx**2 +  # INSERT BC: q[v[i+1,j]]/dx**2
                  (q[v[i,j-1]]-2*q[v[i,j]]+q[v[i,j+1]])/dy**2 
                  )
    # Case 9: interior
    i=np.arange(1,nx-1)
    j=np.arange(1,ny-2)
    #for i in range(1,nx-1):
        #for j in range(1,ny-2):
    ddq[v[1:nx-1,1:ny-2]] = np.add(
                                  (q[v[i-1,1:ny-2]]-2*q[v[i,1:ny-2]]+q[v[i+1,1:ny-2]])/dx**2,
                                  (q[v[1:nx-1,j-1]]-2*q[v[1:nx-1,j]]+q[v[1:nx-1,j+1]])/dy**2 
                                  )  
    
    
    return ddq

#############################################################################################################
def bcu(i,j):
    # NOTE: given a 'u' boundary coordinate, returns scalar value of the boundary condition
    # Currently, this routine enforces no-flow-through at the top,bottom, right, and left walls.
    
    if (hasattr(i, "__len__") and hasattr(j, "__len__")):
        bcon = np.zeros([len(i),len(j)],dtype=float)
    elif hasattr(i, "__len__"):
        bcon = np.zeros([len(i)],dtype=float)
    elif hasattr(j, "__len__"):
        bcon = np.zeros([len(j)],dtype=float)  
    else:
          bcon=0

    return bcon

#############################################################################################################
def bcv(i,j):
    # NOTE: given a 'v' boundary coordinate, returns scalar value of the boundary condition
    # Currently, this routine enforces no-flow-through at the top,bottom, right, and left walls.
    if (hasattr(i, "__len__") and hasattr(j, "__len__")):
        bcon = np.zeros([len(i),len(j)],dtype=float)
    elif hasattr(i, "__len__"):
        bcon = np.zeros([len(i)],dtype=float)
    elif hasattr(j, "__len__"):
        bcon = np.zeros([len(j)],dtype=float)  
    else:
          bcon=0      
    return bcon

#############################################################################################################
def gcu(top,vel):
    # NOTE: given 'u' value of interior velocity, returns ghost cell value
    # The value of the ghost cell is the opposite of the adjacent, interior point for stationary walls.
    # Assuming that the lid moves at uwall, then @ the top boundary v[i+-1,j] = -v[i,j] and u[i,j+-1]=2*uwall-u[i,j]
    # June 4, 2021: should also accept arrays of vel.   
    if top==False:
        gcell=-vel
    elif top==True:
        gcell=2*uwall-vel
    return gcell

#############################################################################################################
def gcv(top,vel):
    # NOTE: given 'v' value of interior velocity, returns ghost cell value
    # The value of the ghost cell is the opposite of the adjacent, interior point for stationary walls.
    # Assuming that the lid moves at uwall, then @ the top boundary v[i+-1,j] = -v[i,j] and u[i,j+-1]=2*uwall-u[i,j]
    # June 4, 2021: should also accept arrays of vel.  
    gcell=-vel
    return gcell

#############################################################################################################
def BC_L(q):
    # Input(q) ~ size velocity
    # Output(ddq) ~ size velocity
    # NOTE: The proper procedure is BC when u[i+-1,j] or v[i,j+-1]. 
    # NOTE: The proper procedure is GC when u[i,j+-1] or v[i+-1,j].       
                  
    ddq = np.zeros([nq],dtype='float')
    
    ## X-DIR
    # Case 1: bottom-left corner
    i=0
    j=0
    ddq[u[i,j]] = (
                  (bcu(i-1,j))/dx**2 +           # INSERT BC: q[u[i-1,j]]/dx**2
                  (gcu(False,q[u[i,j]]))/dy**2   # INSERT GC: q[u[i,j-1]]/dx**2
                  )
    # Case 2: bottom-right corner
    i=nx-2
    j=0
    ddq[u[i,j]] = (
                  (bcu(i+1,j))/dx**2 +           # INSERT BC: q[u[i+1,j]]/dx**2
                  (gcu(False,q[u[i,j]]))/dy**2   # INSERT GC: q[u[i,j-1]]/dy**2
                  )
    
    # Case 5: bottom wall (excluding corners)
    j=0
    i=np.arange(1,nx-2)
    #for i in range(1,nx-2):
    ddq[u[i,j]] = (
                  (gcu(False,q[u[i,j]]))/dy**2   # INSERT GC: q[u[i,j-1]]/dy**2
                  )     
    
    # Case 3: top-left corner ***
    i=0
    j=ny-1
    ddq[u[i,j]] = (
                  (bcu(i-1,j))/dx**2 +            # INSERT BC: q[u[i-1,j]]/dx**2
                  (gcu(True,q[u[i,j]]))/dy**2     # INSERT GC: q[u[i,j+1]]/dy**2
                  )
    # Case 4: top-right corner ***
    i=nx-2
    j=ny-1
    ddq[u[i,j]] = (
                  (bcu(i+1,j))/dx**2 +             # INSERT BC: q[u[i+1,j]]/dx**2
                  (gcu(True,q[u[i,j]]))/dy**2      # INSERT GC: q[u[i,j+1]]/dy**2
                  )
    # Case 6: top wall (excluding corners) ***
    j=ny-1
    i=np.arange(1,nx-2)
    #for i in range(1,nx-2):
    ddq[u[i,j]] = (
                  (gcu(True,q[u[i,j]]))/dy**2   # INSERT GC: q[u[i,j+1]]/dy**2
                  )
    
    # Case 7: left wall (excluding corners)
    i=0
    j=np.arange(1,ny-1)
    #for j in range(1,ny-1):
    ddq[u[i,j]] = (
                  (bcu(i-1,j))/dx**2            # INSERT BC: q[u[i-1,j]]/dx**2
                  )
    
    # Case 8: right wall (excluding corners)
    i=nx-2
    j=np.arange(1,ny-1)
    #for j in range(1,ny-1):
    ddq[u[i,j]] = (
                  (bcu(i+1,j))/dx**2             # INSERT BC: q[u[i+1,j]]/dx**2
                  )
        
    # Case 9: interior (no bc)

    ## Y-DIR       
    # Case 1: bottom-left corner
    i=0
    j=0
    ddq[v[i,j]] = (
                  (gcv(False,q[v[i,j]]))/dx**2 + # INSERT GC: q[v[i-1,j]]/dx**2
                  (bcv(i,j-1))/dy**2             # INSERT BC: q[v[i,j-1]]/dx**2
                  )
    # Case 2: bottom-right corner
    i=nx-1
    j=0
    ddq[v[i,j]] = (
                  (gcv(False,q[v[i,j]]))/dx**2 + # INSERT GC: q[v[i+1,j]]/dx**2
                  (bcv(i,j-1))/dy**2             # INSERT BC: q[v[i,j-1]]/dy**2
                  )
    # Case 5: bottom wall (excluding corners)
    j=0
    i=np.arange(1,nx-1)
    #for i in range(1,nx-1):
    ddq[v[i,j]] = (
                  (bcv(i,j-1))/dy**2         # INSERT BC: q[v[i,j-1]]/dy**2
                  )
        
    # Case 3: top-left corner ***
    i=0
    j=ny-2
    ddq[v[i,j]] = (
                  (gcv(True,q[v[i,j]]))/dx**2 +  # INSERT GC: q[v[i-1,j]]/dx**2
                  (bcv(i,j+1))/dy**2             # INSERT BC: q[v[i,j+1]]/dy**2
                  )
    # Case 4: top-right corner ***
    i=nx-1
    j=ny-2
    ddq[v[i,j]] = (
                  (gcv(True,q[v[i,j]]))/dx**2 +  # INSERT GC: q[v[i+1,j]]/dx**2
                  (bcv(i,j+1))/dy**2             # INSERT BC: q[v[i,j+1]]/dy**2
                  )
    # Case 6: top wall (excluding corners) ***
    j=ny-2
    i=np.arange(1,nx-1)
    #for i in range(1,nx-1):
    ddq[v[i,j]] = (
                  (bcv(i,j+1))/dy**2         # INSERT BC: q[v[i,j+1]]/dy**2
                  ) 
    
    # Case 7: left wall (excluding corners)
    i=0
    j=np.arange(1,ny-2)
    #for j in range(1,ny-2):
    ddq[v[i,j]] = (
                  (gcv(False,q[v[i,j]]))/dx**2   # INSERT GC: q[v[i-1,j]]/dx**2
                  )
    # Case 8: right wall (excluding corners)
    i=nx-1
    j=np.arange(1,ny-2)
    #for j in range(1,ny-2):
    ddq[v[i,j]] = (
                  (gcv(False,q[v[i,j]]))/dx**2   # INSERT GC: q[v[i+1,j]]/dx**2
                  )
    # Case 9: interior (no bc)
    
    return ddq

#############################################################################################################

def Ad(q):
    # Input (q) ~ size velocity
    # Output(n) ~ size velocity
    # boundary condition look-up helper-function: bc
    # NOTE: update the gc at the top boundary (interpolate correctly)
    # Important note: Ad = -Nonlinear
    
    n = np.zeros([nq],dtype='float')
    
    ### Nx ###     
    # Case 1: left wall
    i=0
    j=np.arange(1,ny-1)
    #for j in range(1,ny-1):
    uxl = (bcu(i-1,j)  +q[u[i,j]])/2     # BC: q[u[i-1,j]]
    uxr = (q[u[i,j]]      +q[u[i+1,j]])/2   
    vxt = (q[v[i,j]]      +q[v[i+1,j]])/2
    vxb = (q[v[i,j-1]]    +q[v[i+1,j-1]])/2
    uyt = (q[u[i,j]]      +q[u[i,j+1]])/2
    uyb = (q[u[i,j-1]]    +q[u[i,j]])/2
    dx_uxux = (uxr**2-uxl**2)/dx
    dy_uyvx = (uyt*vxt-uyb*vxb)/dy
    n[u[i,j]] = dx_uxux+dy_uyvx
        
    # Case 2: right wall
    i=nx-2
    j=np.arange(1,ny-1)
    #for j in range(1,ny-1):
    uxl = (q[u[i-1,j]] +q[u[i,j]])/2
    uxr = (q[u[i,j]]   +bcu(i+1,j))/2   # BC: q[u[i+1,j]]
    vxt = (q[v[i,j]]   +q[v[i+1,j]])/2
    vxb = (q[v[i,j-1]] +q[v[i+1,j-1]])/2
    uyt = (q[u[i,j]]   +q[u[i,j+1]])/2
    uyb = (q[u[i,j-1]] +q[u[i,j]])/2
    dx_uxux = (uxr**2-uxl**2)/dx
    dy_uyvx = (uyt*vxt-uyb*vxb)/dy
    n[u[i,j]] = dx_uxux+dy_uyvx
    
    # Case 3: bottom wall
    j=0
    i=np.arange(1,nx-2)
    #for i in range(1,nx-2):
    uxl = (q[u[i-1,j]]    +q[u[i,j]])/2
    uxr = (q[u[i,j]]      +q[u[i+1,j]])/2   
    vxt = (q[v[i,j]]      +q[v[i+1,j]])/2
    vxb = (bcv(i,j-1)+bcv(i+1,j-1))/2          # BC: +q[v[i,j-1]]+q[v[i+1,j-1]]
    uyt = (q[u[i,j]]      +q[u[i,j+1]])/2
    uyb = (gcu(False,q[u[i,j]])+q[u[i,j]])/2      # GC: q[u[i,j-1]]
    dx_uxux = (uxr**2-uxl**2)/dx
    dy_uyvx = (uyt*vxt-uyb*vxb)/dy
    n[u[i,j]] = dx_uxux+dy_uyvx
    
    # Case 4: top wall ***
    j=ny-1
    i=np.arange(1,nx-2)
    #for i in range(1,nx-2):
    uxl = (q[u[i-1,j]]  +q[u[i,j]])/2
    uxr = (q[u[i,j]]    +q[u[i+1,j]])/2   
    vxt = (bcv(i,j)+bcv(i+1,j))/2               # BC: q[v[i,j]]+q[v[i+1,j]]
    vxb = (q[v[i,j-1]]  +q[v[i+1,j-1]])/2
    uyt = (q[u[i,j]]    +gcu(True,q[u[i,j]]))/2 # GC: q[u[i,j+1]]
    uyb = (q[u[i,j-1]]  +q[u[i,j]])/2
    dx_uxux = (uxr**2-uxl**2)/dx
    dy_uyvx = (uyt*vxt-uyb*vxb)/dy
    n[u[i,j]] = dx_uxux+dy_uyvx
    
    # Case 5: bottom-left corner
    i=0
    j=0
    uxl = (bcu(i-1,j)+q[u[i,j]])/2             # BC: q[u[i-1,j]]
    uxr = (q[u[i,j]]      +q[u[i+1,j]])/2   
    vxt = (q[v[i,j]]      +q[v[i+1,j]])/2
    vxb = (bcv(i,j-1)+bcv(i+1,j-1))/2       # BC: q[v[i,j-1]]+q[v[i+1,j-1]]
    uyt = (q[u[i,j]]      +q[u[i,j+1]])/2
    uyb = (gcu(False,q[u[i,j]])+q[u[i,j]])/2   # GC: q[u[i,j-1]]
    dx_uxux = (uxr**2-uxl**2)/dx
    dy_uyvx = (uyt*vxt-uyb*vxb)/dy
    n[u[i,j]] = dx_uxux+dy_uyvx
    
    # Case 6: bottom-right corner
    i=nx-2
    j=0
    uxl = (q[u[i-1,j]]    +q[u[i,j]])/2
    uxr = (q[u[i,j]]      +bcu(i+1,j))/2      # BC: q[u[i+1,j]]
    vxt = (q[v[i,j]]      +q[v[i+1,j]])/2
    vxb = (bcv(i,j-1)+bcv(i+1,j-1))/2      # BC: q[v[i,j-1]]+q[v[i+1,j-1]]
    uyt = (q[u[i,j]]      +q[u[i,j+1]])/2
    uyb = (gcu(False,q[u[i,j]])+q[u[i,j]])/2  # GC: q[u[i,j-1]]
    dx_uxux = (uxr**2-uxl**2)/dx
    dy_uyvx = (uyt*vxt-uyb*vxb)/dy
    n[u[i,j]] = dx_uxux+dy_uyvx
    
    # Case 7: top-left corner ***
    i=0
    j=ny-1
    uxl = (bcu(i-1,j)+q[u[i,j]])/2              # BC: u[i-1,j]
    uxr = (q[u[i,j]]    +q[u[i+1,j]])/2   
    vxt = (bcv(i,j)  +bcv(i+1,j))/2          # BC: q[v[i,j]]+q[v[i+1,j]]
    vxb = (q[v[i,j-1]]  +q[v[i+1,j-1]])/2
    uyt = (q[u[i,j]]    +gcu(True,q[u[i,j]]))/2 # GC: q[u[i,j+1]]
    uyb = (q[u[i,j-1]]  +q[u[i,j]])/2
    dx_uxux = (uxr**2-uxl**2)/dx
    dy_uyvx = (uyt*vxt-uyb*vxb)/dy
    n[u[i,j]] = dx_uxux+dy_uyvx
    
    # Case 8: top-right corner ***
    i=nx-2
    j=ny-1
    uxl = (q[u[i-1,j]]  +q[u[i,j]])/2
    uxr = (q[u[i,j]]    +bcu(i+1,j))/2           # BC: q[u[i+1,j]]
    vxt = (bcv(i,j)  +bcv(i+1,j))/2           # BC: q[v[i,j]]+q[v[i+1,j]]
    vxb = (q[v[i,j-1]]  +q[v[i+1,j-1]])/2
    uyt = (q[u[i,j]]    +gcu(True,q[u[i,j]]))/2  # GC: q[u[i,j+1]]
    uyb = (q[u[i,j-1]]  +q[u[i,j]])/2
    dx_uxux = (uxr**2-uxl**2)/dx
    dy_uyvx = (uyt*vxt-uyb*vxb)/dy
    n[u[i,j]] = dx_uxux+dy_uyvx
    
    if 1: 
        # Case 9: interior
        i=np.arange(1,nx-2)
        j=np.arange(1,ny-1)
        #for i in range(1,nx-2):
            #for j in range(1,ny-1):
        uxl = (q[u[i-1,1:ny-1]] +q[u[i,1:ny-1]])/2
        uxr = (q[u[i,1:ny-1]]   +q[u[i+1,1:ny-1]])/2   
        vxt = (q[v[i,1:ny-1]]   +q[v[i+1,1:ny-1]])/2
        vxb = (q[v[1:nx-2,j-1]] +q[v[1+1:nx-2+1,j-1]])/2
        uyt = (q[u[1:nx-2,j]]   +q[u[1:nx-2,j+1]])/2
        uyb = (q[u[1:nx-2,j-1]] +q[u[1:nx-2,j]])/2
        dx_uxux = (uxr**2-uxl**2)/dx
        dy_uyvx = (uyt*vxt-uyb*vxb)/dy
        n[u[1:nx-2,1:ny-1]] = dx_uxux+dy_uyvx
    
    ### Ny ###    
    # Case 1: left wall
    i=0
    j=np.arange(1,ny-2)
    #for j in range(1,ny-2):
    vyt = (q[v[i,j]]      +q[v[i,j+1]])/2
    vyb = (q[v[i,j-1]]    +q[v[i,j]])/2
    uyl = (bcu(i-1,j)+bcu(i-1,j+1))/2    # BC: q[u[i-1,j+1]]+q[u[i-1,j]]
    uyr = (q[u[i,j]]      +q[u[i,j+1]])/2
    vxl = (gcv(False,q[v[i,j]])+q[v[i,j]])/2   # GC: q[v[i-1,j]]
    vxr = (q[v[i,j]]      +q[v[i+1,j]])/2
    dx_uyvx = (uyr*vxr-uyl*vxl)/dx
    dy_vyvy = (vyt**2-vyb**2)/dy
    n[v[i,j]] =dx_uyvx+dy_vyvy
    
    # Case 2: right wall
    i=nx-1
    j=np.arange(1,ny-2)
    #for j in range(1,ny-2):
    vyt = (q[v[i,j]]    +q[v[i,j+1]])/2
    vyb = (q[v[i,j-1]]  +q[v[i,j]])/2
    uyl = (q[u[i-1,j]]  +q[u[i-1,j+1]])/2
    uyr = (bcu(i,j)  +bcu(i,j+1))/2           # BC: q[u[i,j]]+q[u[i,j+1]]
    vxl = (q[v[i-1,j]]  +q[v[i,j]])/2
    vxr = (q[v[i,j]]    +gcv(False,q[v[i,j]]))/2 # GC: q[v[i+1,j]]
    dx_uyvx = (uyr*vxr-uyl*vxl)/dx
    dy_vyvy = (vyt**2-vyb**2)/dy
    n[v[i,j]] =dx_uyvx+dy_vyvy
    
    # Case 3: bottom wall
    j=0
    i=np.arange(1,nx-1)
    #for i in range(1,nx-1):
    vyt = (q[v[i,j]]      +q[v[i,j+1]])/2
    vyb = (bcv(i,j-1)  +q[v[i,j]])/2       # BC: q[v[i,j-1]]
    uyl = (q[u[i-1,j]]    +q[u[i-1,j+1]])/2
    uyr = (q[u[i,j]]      +q[u[i,j+1]])/2
    vxl = (q[v[i-1,j]]    +q[v[i,j]])/2
    vxr = (q[v[i,j]]      +q[v[i+1,j]])/2
    dx_uyvx = (uyr*vxr-uyl*vxl)/dx
    dy_vyvy = (vyt**2-vyb**2)/dy
    n[v[i,j]] =dx_uyvx+dy_vyvy
    
    # Case 4: top wall
    j=ny-2
    i=np.arange(1,nx-1)
    #for i in range(1,nx-1):
    vyt = (q[v[i,j]]   +bcv(i,j+1))/2   # BC: q[v[i,j+1]]
    vyb = (q[v[i,j-1]] +q[v[i,j]])/2
    uyl = (q[u[i-1,j]] +q[u[i-1,j+1]])/2
    uyr = (q[u[i,j]]   +q[u[i,j+1]])/2
    vxl = (q[v[i-1,j]] +q[v[i,j]])/2
    vxr = (q[v[i,j]]   +q[v[i+1,j]])/2
    dx_uyvx = (uyr*vxr-uyl*vxl)/dx
    dy_vyvy = (vyt**2-vyb**2)/dy
    n[v[i,j]] =dx_uyvx+dy_vyvy
    
    # Case 5: bottom-left corner
    i=0
    j=0
    vyt = (q[v[i,j]]      +q[v[i,j+1]])/2
    vyb = (bcv(i,j-1)  +q[v[i,j]])/2         # BC: q[v[i,j-1]]
    uyl = (bcu(i-1,j)  +bcu(i-1,j+1))/2   # BC: q[u[i-1,j]]+q[u[i-1,j+1]]
    uyr = (q[u[i,j]]      +q[u[i,j+1]])/2
    vxl = (gcv(False,q[v[i,j]])+q[v[i,j]])/2 # GC: q[v[i-1,j]]
    vxr = (q[v[i,j]]      +q[v[i+1,j]])/2
    dx_uyvx = (uyr*vxr-uyl*vxl)/dx
    dy_vyvy = (vyt**2-vyb**2)/dy
    n[v[i,j]] =dx_uyvx+dy_vyvy
    
    # Case 6: bottom-right corner
    i=nx-1
    j=0
    vyt = (q[v[i,j]]      +q[v[i,j+1]])/2
    vyb = (bcv(i,j-1)  +q[v[i,j]])/2         # BC: q[v[i,j-1]]
    uyl = (q[u[i-1,j]]    +q[u[i-1,j+1]])/2
    uyr = (bcu(i,j)    +bcu(i,j+1))/2     # BC: q[u[i,j]]+q[u[i,j+1]]
    vxl = (q[v[i-1,j]]    +q[v[i,j]])/2
    vxr = (q[v[i,j]]+gcv(False,q[v[i,j]]))/2 # GC: q[v[i+1,j]]
    dx_uyvx = (uyr*vxr-uyl*vxl)/dx
    dy_vyvy = (vyt**2-vyb**2)/dy
    n[v[i,j]] =dx_uyvx+dy_vyvy
    
    # Case 7: top-left corner
    i=0
    j=ny-2
    vyt = (q[v[i,j]]      +bcv(i,j+1))/2     # BC: q[v[i,j+1]]
    vyb = (q[v[i,j-1]]    +q[v[i,j]])/2
    uyl = (bcu(i-1,j)  +bcu(i-1,j+1))/2   # BC: q[u[i-1,j]]+q[u[i-1,j+1]]
    uyr = (q[u[i,j]]      +q[u[i,j+1]])/2
    vxl = (gcv(True,q[v[i,j]])+q[v[i,j]])/2  # GC: q[v[i-1,j]]
    vxr = (q[v[i,j]]      +q[v[i+1,j]])/2
    dx_uyvx = (uyr*vxr-uyl*vxl)/dx
    dy_vyvy = (vyt**2-vyb**2)/dy
    n[v[i,j]] =dx_uyvx+dy_vyvy
    
    # Case 8: top-right corner
    i=nx-1
    j=ny-2
    vyt = (q[v[i,j]]    +bcv(i,j+1))/2       # BC: q[v[i,j+1]]
    vyb = (q[v[i,j-1]]  +q[v[i,j]])/2
    uyl = (q[u[i-1,j]]  +q[u[i-1,j+1]])/2
    uyr = (bcu(i,j)  +bcu(i,j+1))/2       # BC: q[u[i,j]]+q[u[i,j+1]]
    vxl = (q[v[i-1,j]]  +q[v[i,j]])/2
    vxr = (q[v[i,j]]+gcv(True,q[v[i,j]]))/2  # GC: q[v[i+1,j]]
    dx_uyvx = (uyr*vxr-uyl*vxl)/dx
    dy_vyvy = (vyt**2-vyb**2)/dy
    n[v[i,j]] =dx_uyvx+dy_vyvy

    if 1: 
        # Case 9: interior
        i=np.arange(1,nx-1)
        j=np.arange(1,ny-2)
        #for i in range(1,nx-1):
            #for j in range(1,ny-2):
        vyt = (q[v[1:nx-1,j]]   +q[v[1:nx-1,j+1]])/2
        vyb = (q[v[1:nx-1,j-1]] +q[v[1:nx-1,j]])/2
        uyl = (q[u[i-1,1:ny-2]] +q[u[i-1,1+1:ny-2+1]])/2
        uyr = (q[u[1:nx-1,j]]   +q[u[1:nx-1,j+1]])/2
        vxl = (q[v[i-1,1:ny-2]] +q[v[i,1:ny-2]])/2
        vxr = (q[v[i,1:ny-2]]   +q[v[i+1,1:ny-2]])/2
        dx_uyvx = (uyr*vxr-uyl*vxl)/dx
        dy_vyvy = (vyt**2-vyb**2)/dy
        n[v[1:nx-1,1:ny-2]] =dx_uyvx+dy_vyvy

    # Recall that the A = - advective term
    n = -n 
    
    return n
 
#############################################################################################################   
def check_matrix():
    
    # Check Divergence Operator (i.e. plot sparsity pattern)
    D = np.zeros([ng,nq])
    for i in range(0,nq):
        z = np.zeros([nq])
        z[i] = 1
        D[:,i] = div(z)
    plt.figure(2)
    plt.spy(D)
    str1 = 'nz=%d'%(np.count_nonzero(D))
    plt.text(60,40,str1)
    ax1 = plt.gca()
    ax1.set_ylabel('ng')
    ax1.set_xlabel('nq') 
    plt.show()
    print('rank div operator =',np.linalg.matrix_rank(D))
    print('nq =',nq)
    print('ng =',ng)
       
    # Check Gradient Operator 
    G = np.zeros([nq,ng])
    for i in range(0,ng):
        z = np.zeros([ng])
        z[i] = 1
        G[:,i] = grad(z)
    
    plt.figure(3)
    plt.spy(G)
    str1 = 'nz=%d'%(np.count_nonzero(G))
    plt.text(60,40,str1)
    ax1 = plt.gca()
    ax1.set_ylabel('nq')
    ax1.set_xlabel('ng') 
    plt.show()
    print('rank grad operator =',np.linalg.matrix_rank(G))
    print('nq =',nq)
    print('ng =',ng)

    # By D = -G^T, the number of nonzero elements should be the same.
    # Therefore, the inifinity norm of D+(-G) should be zero. 
    print('D+G^T=0? D+G^T=%d'%(np.linalg.norm(D+G.T,np.inf)))
    # confirmed 05/02/2021
    
    # Check Laplacian Operator
    L = np.zeros([nq,nq])
    for i in range(0,nq):
        z = np.zeros([nq])
        z[i] = 1
        L[:,i] = lap(z)
    plt.figure(4)
    plt.spy(L)
    str1 = 'nz=%d'%(np.count_nonzero(L))
    plt.text(60,40,str1)
    ax1 = plt.gca()
    ax1.set_ylabel('nq')
    ax1.set_xlabel('nq') 
    plt.show()
    
    # By L = L^T, the infinity norm of L-L^T should be zero. 
    print('L=L^T? L-L^T=%d'%(np.linalg.norm(L-L.T,np.inf)),'\n')
    # confirmed 05/02/2021
    
    # Check Nonlinear/Advective Operator (A)
    N = np.zeros([nq,nq])
    for i in range(0,nq):
        z = np.zeros([nq])
        z[i] = 1
        N[:,i] = Ad(z)
    plt.figure(5)
    plt.spy(N)
    str1 = 'nz=%d'%(np.count_nonzero(N))
    plt.text(60,40,str1)
    ax1 = plt.gca()
    ax1.set_ylabel('nq')
    ax1.set_xlabel('nq')
    plt.show()
    
    # Optional: Check Positive-Definite Mom-Eq
    if 0:
        M = np.zeros([nq,nq])
        for i in range(0,nq):
            z = np.zeros([nq])
            z[i] = 1
            M[:,i] = R(z)
        print('mom-eq eigvals =',np.linalg.eigvals(M))
        print('mom-operator symmetric =',check_symmetric(M))
        
    # Optional: Check Positive-Definite and Symmetric PP-Eq
    if 0:
        M = np.zeros([ng,ng])
        i = 0
        z = np.zeros([ng])
        M[:,i] = -div(Rinv(grad(z)))
        for i in range(1,ng):
            z = np.zeros([ng])
            z[i] = 1
            M[:,i] = -div(Rinv(grad(z)))
        print('pp-eq eigvals =',np.linalg.eigvals(M))
        print('PP-operator symmetric =',check_symmetric(M))
    if 0:
        M = np.zeros([ng,ng])
        i = 0
        z = np.zeros([ng])
        M[:,i] = -100
        for i in range(1,ng):
            z = np.zeros([ng])
            z[i] = 1
            M[:,i] = -div(grad(z))
        print('pp-eq eigvals -div(grad(z)) =',np.linalg.eigvals(M))
        print('PP-operator symmetric =',check_symmetric(M))
        
#############################################################################################################  
def check_symmetric(a, rtol=1e-05, atol=1e-08):
    return np.allclose(a, a.T, rtol=rtol, atol=atol)
        
#############################################################################################################  
#############################################################################################################  
#############################################################################################################  
import numpy as np
import matplotlib.pyplot as plt
import time
import pickle
import os
from datetime import datetime

NS_solver()
    
