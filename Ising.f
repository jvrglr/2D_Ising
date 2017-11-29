      program ISING_2D_metropolis
        !position (x,y) in grid is labbeled with with integer i=(y-1)*L+x
        parameter (L=6, N=L*L)
        integer*4 neigh(4,N),s(N)
        call dran_ini(1994)
        open(unit=1, file="data.dat", status="unknown")
        neigh=0
        call Create_neigbours(neigh,L,N)

        do i=1,N
          write(1,*) neigh(2,i)
        enddo

        close(1)
      end

      subroutine Create_neigbours(neigh,L,N)
        !neigh(i,j)= neighbour i of site j
        ! i=1--> right neighbour
        ! i=2--> left neighbour
        ! i=3--> down neighbour
        ! i=4--> up neighbour
        !The implementations of the (periodic) Boundary Conditions is included in this array
        integer*4 neigh(4,N)
        integer*4 k,m

        do k=1,L !coordinate x
          do m=1,L !coordinate y
            j=(m-1)*L+k

            ! right neighbours
            neigh(1,j)=j+1
            if (k.eq.L) then
              neigh(1,j)=(m-1)*L+1 !first coulmn
            endif

            !left neighbours
            neigh(2,j)=j-1
            if (k.eq.1) then
              neigh(2,j)=(m-1)*L+L !last coulmn
            endif

            !down neighbours
            neigh(3,j)=j-L
            if (m.eq.1) then
              neigh(3,j)=k+L*(L-1) !last row
            endif

            !up neighbpurs
            neigh(4,j)=j+L
            if (m.eq.L) then
              neigh(4,j)=k !first row BC
            endif

          enddo
        enddo

        return
      end subroutine
