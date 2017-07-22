import numpy as np
cimport numpy as np
cimport cython

cdef extern from "time.h" nogil:
    ctypedef int time_t
    time_t time(time_t*)

def computeQe(dict slicelist, list node_data, int discrsize, float drate, float trate, float lrate, float stemlen):
    return c_computeQe(slicelist, node_data, discrsize, drate, trate, lrate, stemlen)

def nodeLimitter(genetree, int discrsize, int leafslice):
    return c_nodeLimitter(genetree, discrsize, leafslice)

def computeProb(rateDens, genetree not None, dict name2node, dict rankedge, list node_data, int discrsize, float stemlen, float drate, float trate, np.ndarray[np.float_t, ndim=4] Qef):
    return c_computeProb(rateDens, genetree, name2node, rankedge, node_data, discrsize, stemlen, drate, trate, Qef)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.ndarray[np.float_t, ndim=4] c_computeQe(dict slicelist, list node_data, int discrsize, float drate, float trate, float lrate, float stemlen):
    # should return a numpy array
    cdef:
        int ctmp_size
        int t
        int s
        float h
        int e_ind
        int f_ind
        float wnorm
        int n_edges = len(node_data)
        float qehf_sum
        float rates = drate + trate + lrate
    # link view to numpy array
    # Q_e should contains the value for t_end (parent node for the edge e)
    cdef np.ndarray[np.float_t, ndim=2, mode='c'] Qe = np.zeros((n_edges, discrsize+1), dtype=np.float)
    cdef np.ndarray[np.float_t, ndim=4, mode='c'] Qef = np.zeros((n_edges, n_edges, discrsize+1, discrsize+1), dtype=np.float)

    cdef np.ndarray[np.float_t, ndim=2, mode='c'] Qe_t, Qe_k1, Qe_k2, Qe_k3, Qe_k4, Qef_k1, Qef_k2, Qef_k3, Qef_k4
    cdef np.ndarray[np.float_t, ndim=3, mode='c'] Qef_st
    cdef np.ndarray[np.float_t, ndim=2, mode='c'] tmpmask
    
    for rank in xrange(max(slicelist.keys()), -1, -1):
        edgelist = slicelist[rank]
        # corresponding node is node_data[edgelist[0]]
        
        tmpnode = node_data[edgelist[0]]
        t_start = tmpnode.time        
        if tmpnode.up:
            t_end = tmpnode.up.time
        else:
            t_end = t_start + stemlen
        
        ctmp_size = len(edgelist)
        wnorm = 0.0
        if ctmp_size -1 > 0:
            wnorm = trate*1.0 / (ctmp_size-1)

        h = (t_end - t_start)*1.0/discrsize # important since time step
        # can change between two different slice

        Qe_t = np.zeros((ctmp_size, discrsize+1), dtype=np.float)
        Qe_k1 = np.zeros((ctmp_size, discrsize), dtype=np.float)
        Qe_k2 = np.zeros((ctmp_size, discrsize), dtype=np.float)
        Qe_k3 = np.zeros((ctmp_size, discrsize), dtype=np.float)
        Qe_k4 = np.zeros((ctmp_size, discrsize), dtype=np.float)
        
        Qef_st = np.zeros((ctmp_size**2, discrsize+1, discrsize+1), dtype=np.float)
        
        tmpmask = np.ones((ctmp_size,ctmp_size))
        np.fill_diagonal(tmpmask,0)
        # initial condition
        for e_ind, edge_e in enumerate(edgelist):           
            node_e = node_data[edge_e]
            # basic initial condition
            if node_e.is_leaf():
                # this case is expected only when t_start = 0
                assert t_start == 0, "t_start should be 0 at the starting"
                Qe_t[e_ind, 0] = 0
            elif len(node_e.get_children()) == 1:
                Qe_t[e_ind, 0] =  Qe[node_e.children[0].edge_i, discrsize]
            else:
                Qe_t[e_ind, 0] = Qe[node_e.children[0].edge_i, discrsize] *  Qe[node_e.children[1].edge_i, discrsize]  
            
            # initialisation of Qef(t,t) if e==f
            for t in xrange(discrsize+1):
                Qef_st[e_ind*ctmp_size + e_ind, t, t] = 1

        # discretisation to compute Qe(t)
        for t in xrange(discrsize):
            # t start at 0 which was already computed
            rowsum = np.sum(Qe_t[:, t])
            # Qe_t was computed for all edges under the initial conditions
            # here we compute Qeki for all edges
            qe_t = Qe_t[:, t]
            # compute k1
            Qe_k1[:, t] = drate*(qe_t**2) + lrate - rates*qe_t + wnorm*(rowsum - qe_t)*qe_t 
            # compute k2
            qe_t = Qe_t[:, t] + 0.5*h*Qe_k1[:, t]
            Qe_k2[:, t] = drate*(qe_t**2) + lrate - rates*qe_t + wnorm*(rowsum - qe_t)*qe_t 
            # compute k3
            qe_t = Qe_t[:, t] + 0.5*h*Qe_k2[:, t]
            Qe_k3[:, t] = drate*(qe_t**2) + lrate - rates*qe_t + wnorm*(rowsum - qe_t)*qe_t 
            # compute k4
            qe_t = Qe_t[:, t] + h*Qe_k3[:, t]
            Qe_k4[:, t] = drate*(qe_t**2) + lrate - rates*qe_t + wnorm*(rowsum - qe_t)*qe_t
            # compute Qe_t        
            Qe_t[:,t+1] = Qe_t[:, t] + h*(Qe_k1[:, t] + 2*Qe_k2[:, t] + 2*Qe_k3[:, t] + Qe_k4[:, t])/6.0

        # this value will change at each loop
        # so no need to initialize multiple time
        Qef_k1 = np.zeros((ctmp_size**2, discrsize), dtype=np.float)
        Qef_k2 = np.zeros((ctmp_size**2, discrsize), dtype=np.float)
        Qef_k3 = np.zeros((ctmp_size**2, discrsize), dtype=np.float)
        Qef_k4 = np.zeros((ctmp_size**2, discrsize), dtype=np.float)
        for t in xrange(discrsize):
            # for s > t
            for s in xrange(t, discrsize):
                qe_s = Qe_t[:, s] # this is the value for Qe(s), which is Qe(t) at the start
                qef_st = Qef_st[:, t, s] # this is the previous value
                # of Qef(s,t)
                # for the first value, s=t and Qef(s,t) = {0,1}
                qef_st_r = qef_st.reshape((ctmp_size, ctmp_size))
                qeqef_st = ((qe_s*qef_st_r).T).ravel()
                
                # compute k1
                Qef_k1[:, s] = 2*drate*qeqef_st - rates*qef_st + wnorm*((np.dot(tmpmask, qef_st_r)*qe_s[:, None]).ravel() + qef_st*sum(qe_s) - qeqef_st)
                
                #  compute k2
                qef_st = Qef_st[:, t, s] + 0.5*h*Qef_k1[:, s]               
                Qef_k2[:, s] = 2*drate*qeqef_st - rates*qef_st + wnorm*((np.dot(tmpmask, qef_st_r)*qe_s[:, None]).ravel() + qef_st*sum(qe_s) - qeqef_st)
                
                # compute k3
                qef_s = Qef_st[:, t, s] + 0.5*h*Qef_k2[:, s]
                Qef_k3[:, s] = 2*drate*qeqef_st - rates*qef_st + wnorm*((np.dot(tmpmask, qef_st_r)*qe_s[:, None]).ravel() + qef_st*sum(qe_s) - qeqef_st)
                
                # compute k4
                qef_s = Qef_st[:, t, s] + h*Qef_k2[:, s]
                Qef_k4[:, s] = 2*drate*qeqef_st - rates*qef_st + wnorm*((np.dot(tmpmask, qef_st_r)*qe_s[:, None]).ravel() + qef_st*sum(qe_s) - qeqef_st)
            
                Qef_st[:,t, s+1] = Qef_st[:, t, s] + h*(Qef_k1[:, s] + 2*Qef_k2[:, s] + 2*Qef_k3[:, s] + Qef_k4[:, s])/6.0

        for e_ind, edge_e in enumerate(edgelist): 
            # this is the value at t_end (parent node) of the edge
            # for leaves, Qe[edge] is expected to be 0 
            #print Qe[edge_e, :].size
            #print Qe_t[e_ind, :].size
            Qe[edge_e, :] = Qe_t[e_ind, :]
            for f_ind, edge_f in enumerate(edgelist): 
                Qef[edge_e, edge_f, :, :] = Qef_st[e_ind*ctmp_size + f_ind, :, :]
            
        # update Qef for  --e-->v--f-->
        # start by finding the speciation time
        snode, g_edge = None, None # for time slice
        all_leaves = True
        for edge in edgelist:
            corr_node = node_data[edge]
            if len(corr_node.get_children()) == 2:
                snode, g_edge = corr_node, edge # only node with two children so just use first
                break
            all_leaves = (all_leaves and corr_node.is_leaf())
        # speciation node found
        # meaning this is not the leaves time slice
        # if the speciation node g is not found
        # then something is wrong with the slice
        if snode:
            gp_edge, gpp_edge = [child.edge_i for child in snode.get_children()]
            for e_edge in edgelist:
                for below_rank in xrange(max(slicelist.keys()), rank, -1):
                    # f is the direct child of e
                    for f_edge in slicelist[below_rank]:
                        # outgoing edge interval discr (f)
                        # speciation time is either 0 or 10. 
                        for t in xrange(discrsize+1):
                            # incoming edge (e) discr
                            for s in xrange(discrsize+1):
                                # sum over contemporain to g edge
                                qehf_sum = sum([Qef[e_edge, int(h_edge), 0, s]*Qef[int(h_edge), f_edge, t, discrsize] for h_edge in edgelist if h_edge != g_edge])
                                # compute rest
                                Qef[e_edge, f_edge, t, s] =  qehf_sum + Qef[e_edge, g_edge, 0, s]*(Qef[gp_edge, f_edge, t, discrsize]*Qe[gpp_edge, discrsize] + Qef[gpp_edge, f_edge, t, discrsize]*Qe[gp_edge, discrsize])
        elif not all_leaves:
            #print [node_data[edg] for edg in edgelist]
            raise ValueError("time slice without a speciation delimiting it")
    return Qef


cpdef int get_discr_size(int gsize, int discrsize, int leafslice):
    while (discrsize-1)*leafslice < gsize:
        discrsize += 1
    return discrsize

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void c_nodeLimitter(genetree, int discrsize, int leafslice):
    # discrsize should be at list Number_of_internal + 1
    # this will attempt to find the genetree to species tree node mapping
    # bounds based on topology only
    cdef: 
        int lowSlice
        int lowTime

    for post, node in genetree.iter_prepostorder():
        if node.is_leaf():
            # leaves can only be mapped to first discr
            node.add_features(lowTime=0) # starting time at bottom
            node.add_features(upTime=0) # starting time at top
            node.add_features(lowSlice=leafslice)
            node.add_features(upSlice=leafslice)
        else:
            if not post and node.is_root():
                node.add_features(upSlice=0)
                # set upper limit to just before the specietree tree stem tip
                node.add_features(upTime=discrsize-1)
                
            elif post:
                child1, child2 =  node.get_children()
                lowSlice = min(child1.lowSlice, child2.lowSlice)
                lowTime = max(child1.lowTime , child2.lowTime) + 1
                if lowTime >= discrsize:
                    lowTime -= discrsize
                    lowSlice -= 1
                node.add_features(lowSlice=lowSlice)
                node.add_features(lowTime=lowTime)
            
            else:
                upSlice = node.up.upSlice
                upTime = node.up.upTime - 1
                if upTime < 0 : 
                    upTime = discrsize-1
                    upSlice += 1
                node.add_features(upSlice=upSlice)
                node.add_features(upTime=upTime) 


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef np.float_t [:, :, :] c_computeProb(rateDens, genetree, dict name2node, dict rankedge, list node_data, int discrsize, float stemlen, float drate, float trate, np.ndarray[np.float_t, ndim=4] Qef):
                        
    cdef int gsize = len([gnode for gnode in genetree.traverse()])
    cdef int tot_edge = len(node_data)
    cdef int maxrank = int(max(rankedge.keys()))    

    cprobAx = np.zeros((gsize, tot_edge, discrsize), dtype=np.float)
    cprobSe = np.zeros((gsize, tot_edge, discrsize), dtype=np.float)
    #cdef np.float_t cprobAx[gsize][tot_edge][discrsize]
    #cdef double cprobSe[gsize][tot_edge][discrsize]
    cdef np.float_t [:, :, :] prob_Ax = cprobAx
    cdef np.float_t [:, :, :] prob_Se = cprobSe
    # start by labelling each node of the genetree
    # set species should have been
    # genetree should be labelled also
    cdef:
        int upLim
        float sigma
        float t_x
        float t_z
        float dt
        int lowLim
        float transtmp
        float duptmp
        int e
        int e_discr
        int linked_s_node
        int t0
        int t1
        int t00
        int t01
        int t02
        # list genelist = list(reversed(list(genetree.traverse("levelorder"))))

    for gnode in genetree.traverse("postorder"):

        # we can use this occassion to compute s(e, u)
        gchild = [ch.ind for ch in gnode.get_children()]
        # print 'succcesss 1', gnode.is_leaf()
        if gnode.is_leaf():
            t0 = time(NULL)
            linked_s_node = name2node[gnode.species]

            # gnode.dist is distance to parent
            # get the species node linked to gnode
            # here we are trying to compute s(e,u) = p11(e, z)*sigma
            # for each edge e
            # remember that only the node below is use for an edge
            for edgelist in [rankedge[x] for x in sorted(rankedge.keys(), reverse=True)]:
                #for e in rankedge[gnode.upSlice]:
                
                for e in edgelist:
                    if node_data[e].is_root():
                        dt = stemlen*1.0
                    else:
                        dt = (node_data[e].up.time - node_data[e].time)*1.0                
                    
                    for e_discr in range(discrsize):
                        # go with one round and compute for all leaves before
                        # e = <x, y> ==> t(x) = first discr not at zero
                        #t_x = get_time(e, e_discr, node_data)
                        t_x = node_data[e].time + (e_discr+1)*dt/discrsize
                        sigma = rateDens.pdf(gnode.dist/t_x)
                        # child = [ch.edge_i for ch in node_data[e].get_children()]
                        # for a leaf lowTime and upTime are the same
                        # obviously, Qef is 0 for node that are neither contemporain
                        # nor adjacent. 
                        # e = <y, z>
                        # linked_s_node = x and f = <xx, x>
                        prob_Se[gnode.ind, e, e_discr] = Qef[e, linked_s_node, 0, e_discr+1]*sigma
                        #print 'mapping node'
                        #print linked_s_node
                        #print 'gnode'
                        #print gnode
                        #print 'snode'
                        #print node_data[e], "\t\t", e_discr
                        #print prob_Se[gnode.ind, e, e_discr]
                        # compute prob_Ax ==> never executed
            t1 = time(NULL)
            print(gnode)
            print("finished in %f"%((t1 - t0)))
        else:
            # for gnode ==> use up time and low time
            # start by finding the time slice and the epoch
            gls = gnode.lowSlice
            glt = gnode.lowTime
            # print '#################'
            # print gnode.upSlice, gnode.upTime
            # print gnode.lowSlice, gnode.lowTime
            # print '##################'
            t0 = time(NULL)

            while gls > gnode.upSlice or (gls==gnode.upSlice and glt <= gnode.upTime):
                print '**********'

                print(glt, gls)
                print gnode
                t00 = time(NULL)
                for e in rankedge[gls]:
                    transtmp = 0.0
                    duptmp = 0.0
                    child = [ch.edge_i for ch in node_data[e].get_children()]
                    # compute prob_Ax
                    if glt == 0:
                        # edge is contempory to a speciation:
                        if len(node_data[e].get_children()) > 1:
                            # edge is the  speciation
                            #if prob_Ax[gnode.ind, e, glt]!=0.0:
                            prob_Ax[gnode.ind, e, glt] = prob_Se[gchild[0], child[0], discrsize-1]*prob_Se[gchild[1], child[1], discrsize-1] + prob_Se[gchild[0], child[1], discrsize-1]*prob_Se[gchild[1], child[0], discrsize-1]
                    else:
                        # edge is not contempory to any speciation
                        if len(rankedge[gls]) > 1:
                            for ctmp_edge in rankedge[gls]:
                                if ctmp_edge != e:
                                    transtmp += prob_Se[gchild[0], e, glt-1]*prob_Se[gchild[1], ctmp_edge, glt-1] + prob_Se[gchild[1], e, glt-1]*prob_Se[gchild[0], ctmp_edge, glt-1]
                            transtmp *= trate*1.0 / (len(rankedge[gls]) -1)
                        duptmp = prob_Se[gchild[0], e, glt-1]*prob_Se[gchild[1], e, glt-1]
                        prob_Ax[gnode.ind, e, glt] = 2*drate*duptmp + transtmp
                t01 = time(NULL)
                print "first loop took: ", (t01-t00)
                for e in rankedge[gls]:
                    transtmp = 0.0
                    duptmp = 0.0
                    child = [ch.edge_i for ch in node_data[e].get_children()]
                    if node_data[e].is_root():
                        dt = stemlen*1.0
                    else:
                        dt = (node_data[e].up.time - node_data[e].time)*1.0

                    t_x = node_data[e].time + (glt+1)*dt/discrsize
                    # find list of z 
                    zlist = []
                    zlt = glt
                    zls = gls
                    while zls <= maxrank and zlt >= 0:
                        # glt is included, since edge represente < glt, glt+1 >
                        if zls < maxrank or (zls==maxrank and zlt>=0):
                            zlist.extend([(zedge,zls,zlt) for zedge in rankedge[zls]])
                        zlt -= 1
                        if zlt < 0:
                            zls += 1
                            zlt += discrsize
                        
                    for zedge, zls, zlt in reversed(zlist):
                        # this start to look repetitif 
                        if node_data[zedge].is_root():
                            dt = stemlen*1.0
                        else:
                            dt = (node_data[zedge].up.time - node_data[zedge].time)*1.0
                        t_z = node_data[zedge].time + zlt*dt/discrsize
                        sigma = rateDens.pdf(gnode.dist/(t_x - t_z))
                        # p11[e, z] = Qef[e, zedge, t_z, t_x] # i correct
                        # e = <x, y> , t_x is glt and t_y is glt+1, assuming that glt max value is discrsize -1 
                        #print '*** Qef', Qef[e, zedge, zlt, glt+1]
                        prob_Se[gnode.ind, e, glt] += Qef[e, zedge, zlt, glt+1]*sigma*prob_Ax[gnode.ind, zedge, zlt]
                        
                        #print '*** a(x,u)',zedge, zls, zlt, prob_Ax[gnode.ind, zedge, zlt], Qef[e, zedge, zlt, glt+1]
                        #print node_data[zedge]
                t02 = time(NULL)
                print "second loop took: ", (t02-t01)
                glt += 1
                if glt >= discrsize:
                    gls -= 1
                    glt -= discrsize        
            t1 = time(NULL)
            print "Node done in ", (t1- t0)
    return prob_Ax
