      program ISING_2D_metropolis
        parameter (L=32, N=L*L)
        integer*4 neigh(4,N),s(N)
        call dran_ini(1994)
        
        write(*,*) L,N
      end
