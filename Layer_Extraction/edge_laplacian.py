'''
NOTE: This code came from recovery.py, which can be found on GitHub:
        https://github.com/yig/harmonic_interpolation
'''

from numpy import *

def gen_edge_laplacian(rows, cols, edge ):
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

    def sum_ij_edge(i,j):
        sum=0
        for h in range(-1,2):
            for w in range(-1,2):
                if edge[i+h][j+w]<255 :
                    sum+=1.0
        return sum-1.0

    Adj = []
    AdjValues = []
    

    for i in xrange( 1, rows-1 ):
        for j in xrange( 1, cols-1 ):

            if edge[i][j]<255:
                nebor_num=sum_ij_edge(i,j)

                if nebor_num>0:
                    assign_value=1/nebor_num
                    ind00 = ij2ind( i,j )

                    count=0
                    for h in range(-1,2):
                        for w in range(-1,2):
                            if (h != 0 or w !=0) and edge[i+h][j+w]<255 :
                                count =count+1
                                indp0 = ij2ind( i+h,j+w )
                                Adj.append( ( ind00, indp0 ) )
                                AdjValues.append(assign_value )






    ## Build the adjacency matrix.
    AdjMatrix = sparse.coo_matrix( ( AdjValues, asarray( Adj ).T ), shape = ( rows*cols, rows*cols ) )
    Mass = sparse.coo_matrix((asarray(AdjMatrix.sum(1)).ravel(), (range(rows * cols), range(rows * cols))))
    L = (Mass - AdjMatrix)
    ## The rows should sum to 0.
    assert (abs(asarray(L.sum(1)).ravel()) < 1e-5).all()


    return L.tocsr()
