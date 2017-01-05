subroutine distance_NN(train,test,n,m,test_row,k,class_NN,trainC,rank,p,chunk)
   use mpi
   implicit none
   integer                :: i, j, k, n, m,copy, train_row, test_row, col, chunk
   real,intent(IN)        :: train(n,m), test(n,m)
   integer,intent(IN)     :: trainC(n)
   integer,intent(out)    :: class_NN(k)
   real                   :: dist(k),d 
   integer                :: p,rank,ierror,status(mpi_status_size),temp2,dummy(n)
   real                   :: gatherelements(n),temp
   real                   :: recvbuf(chunk,m),ppdist(chunk),b(chunk)
   integer,parameter      :: rl = mpi_real,comm = mpi_comm_world

   
!initialise the arrays
if(rank==0)then
   do i  = 1, k
      dist(i)     = 1.E+15
      class_NN(i) = 0
   enddo
endif
call mpi_scatter(train,chunk*m,mpi_real,recvbuf,chunk*m,mpi_real,0, comm ,ierror)

 do train_row = 1,chunk
      d = 0.0
      do col = 1, m
         d = d + (recvbuf(train_row,col) - test(test_row,col))**2
    enddo
     ppdist(train_row)=d
 enddo

call mpi_gather(ppdist(1),chunk,mpi_real,gatherelements(1),chunk, mpi_real,0,comm ,ierror)

if(rank==0)then

 do i=1,n
  dummy(i)=trainC(i)
 enddo

!SORTING
 do i =1,n-1
    do j =i+1,n
       if(gatherelements(i)>gatherelements(j))then
       temp              = gatherelements(j)
       temp2             = dummy(j)
       gatherelements(j) = gatherelements(i)
       dummy(j)          = dummy(i)
       gatherelements(i) = temp
       dummy(i)         = temp2
       endif
    enddo
 enddo

!print*,dummy
 do copy=1,k
  class_NN(copy)  = dummy(copy)
 enddo

endif

end subroutine distance_NN
