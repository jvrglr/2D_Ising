      program ISING_2D_metropolis
        !Simulation of 2D Ising model with zero external field using Metropolis algorithm
        !Bibliography:
        !**************************************************************
        !Toral, R., & Colet, P. (2014). Stochastic numerical methods:
        !an introduction for students and scientists. John Wiley & Sons.
        !Chapter 5
        !
        !Tomás Sintes's notes:
        !https://ifisc.uib-csic.es/users/tomas/MFCS/SquareLatticeIsing.pdf
        !**************************************************************

        !position (x,y) in grid is labbeled with with integer i=(y-1)*L+x
        implicit none
        integer*4 L,N
        parameter (L=100, N=L*L)
        double precision T,u,mag
        double precision hamil,dran_u !Real functions
        double precision h(-4:4),E
        integer*4 i_dran,sum_n !Integer function
        integer*4 neigh(4,N),s(N)
        integer*4 i,j,k,m,rd,steps,step_therm,prod_spin,aux
        call dran_ini(1994)
        open(unit=1, file="averages.dat", status="unknown")
        neigh=0
        call Create_neigbours(neigh,L,N)
        call random_IC(s,N) !random initial condition

        do k=1,10
        T=5.0d0-0.49d0*dble(k) !temperature in J/Kb units

          !Compute possible acceptance h(i)=min(1,exp(-beta*DH({s})))
          !As there are only 4 NN and spin values are -1 or +1--> 5 different values of h
          !h(-4),h(-2),h(0),h(2),h(4)
          do i=-4,0,2
            h(i)=1
          enddo
          h(2)=dexp(-4.0d0/T)
          h(4)=dexp(-8.0d0/T)

          !THERMALIZATION
          steps=1000*N !1000 MC steps
          do i=1,steps
            call rejection(N,s,neigh,h)
          enddo

          !INITIALIZE
          mag=0.0d0

          !MAIN PROGRAM
          !+++++++++++++++++++++++++++++++++++++++++++++++++++
          steps=1000 !Number of measures
          do i=1,steps
            ! This loop pretends to decrease the correlation between measures
            do j=1,N  !ONE mc step
              call rejection(N,s,neigh,h)
            enddo

            !MEASURE
            mag=mag+abs(dble(sum(s)))
          enddo
          mag=mag/(dble(steps)*dble(N)) !So magnetization is in [0,1]
          write(1,*) T,mag

        enddo
        close(1)
      end

      subroutine rejection(N,s,neigh,h)
        integer*4 N,rd,prod_spin,sum_n,i_dran
        integer*4 s(N),neigh(4,N)
        double precision dran_u,u,h(-4:4)
        !MAKES REJECTION STEP
        rd=i_dran(N) !choose spin at random
        prod_spin=s(rd)*sum_n(rd,N,s,neigh)
        u=dran_u()
        if (h(prod_spin).ge.u) then !rejection step (CAN BE IMPROVED)
          s(rd)=-s(rd)
        endif
        return
      end

      integer*4 function sum_n(rd,N,s,neigh)
        integer*4 rd,N,i
        integer*4 s(N),neigh(4,N)
        !sum spins of all nearest neighbours of site rd
        sum_n=0
        do i=1,4
          sum_n=sum_n+s(neigh(i,rd))
        enddo
        return
      end

      double precision function hamil(s,N,T,neigh)
        !compute Hamiltonian.
        !Actually, I'm computing -Hamiltonian/J
        double precision T
        integer*4 N
        integer*4 s(N),neigh(4,N)
        integer*4 i,j
        hamil=0.0d0
        do i=1,N
          do j=1,4 !4=number of NN
            hamil=hamil+s(i)*s(neigh(j,i))
          enddo
        enddo
        return
      end function

      subroutine random_IC(s,N)
        integer*4 i,N,rd
        integer*4 s(N)
        !RANDOM INITIAL CONDITION
        do i=1,N
          rd=i_dran(2)
          s(i)=2*rd-3 !random number, 1 or -1
        enddo
        return
      end

      subroutine Create_neigbours(neigh,L,N)
        !neigh(i,j)= neighbour i of site j
        ! i=1--> right neighbour
        ! i=2--> left neighbour
        ! i=3--> down neighbour
        ! i=4--> up neighbour
        !The implementations of the (periodic) Boundary Conditions is included in this array
        integer*4 k,m,N,L
        integer*4 neigh(4,N)


        do k=1,L !coordinate x (coulmns)
          do m=1,L !coordinate y (rows)
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
