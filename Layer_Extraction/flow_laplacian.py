'''
NOTE: This code came from recovery.py, which can be found on GitHub:
        https://github.com/yig/harmonic_interpolation
'''

from numpy import *

def gen_flow_laplacian(grad,theta, rows, cols, cut_edges = None ):
    '''
    The same as 'gen_symmetric_grid_laplacian1()', except boundary weights are correct.
    
    tested
    (see also test_cut_edges())
    
    '''

    import  numpy
    assert rows > 0
    assert cols > 0

    from scipy import sparse
    
    N = rows
    M = cols
    def ind2ij( ind ):
        assert ind >= 0 and ind < N*M
        return ind // M, ind % M
    def ij2ind( i,j ):
        assert i >= 0 and i < N and j >= 0 and j < M
        return i*M + j

    def ij_flow(i,j):
        assert i >= 0 and i < N and j >= 0 and j < M
        sobel_x= numpy.array([[-1, 0, 1],[-2, 0, 2],[-1, 0, 1]])
        sobel_y= numpy.array([[-1, -2, -1],[0, 0, 0],[1, 2, 1]])
        flow=cos(theta[i][j])*sobel_x-sin(theta[i][j])*sobel_y
        #grad[i][j]=grad[i][j]+0.01
        return flow*grad[i][j]

    Adj = []
    AdjValues = []
    

    for i in xrange( 1, rows-1 ):
        for j in xrange( 1, cols-1 ):

            flow=ij_flow(i,j)

            ind00 = ij2ind( i,j )
            indp0 = ij2ind( i-1,j-1 )
            Adj.append( ( ind00, indp0 ) )
            AdjValues.append( flow[0][0] )

            indp1 = ij2ind(i-1, j)
            Adj.append((ind00, indp1))
            AdjValues.append(flow[0][1])

            indp2 = ij2ind(i-1, j+1)
            Adj.append((ind00, indp2))
            AdjValues.append(flow[0][2])

            indp3 = ij2ind(i, j-1)
            Adj.append((ind00, indp3))
            AdjValues.append(flow[1][0])

            indp4 = ij2ind(i, j)
            Adj.append((ind00, indp4))
            AdjValues.append(flow[1][1])

            indp5 = ij2ind(i, j+1)
            Adj.append((ind00, indp5))
            AdjValues.append(flow[1][2])


            indp6 = ij2ind(i + 1, j-1)
            Adj.append((ind00, indp6))
            AdjValues.append(flow[2][0])

            indp7 = ij2ind(i + 1, j)
            Adj.append((ind00, indp7))
            AdjValues.append(flow[2][1])

            indp8 = ij2ind(i + 1, j+1)
            Adj.append((ind00, indp8))
            AdjValues.append(flow[2][2])

    
    ## Build the adjacency matrix.
    AdjMatrix = sparse.coo_matrix( ( AdjValues, asarray( Adj ).T ), shape = ( rows*cols, rows*cols ) )

    L = AdjMatrix

    return L.tocsr()
