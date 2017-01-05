subroutine nearest_neighbours(train,test,n,m,k,class_assgn,trainC,p1,p2,rank,p)
   use mpi
   implicit none
   integer              ::  test_row, i, k, n, m, count_1, count_2
   real   , intent(in)  ::  train(n,m), test(n,m)
   integer, intent(in)  ::  trainC(n)
   integer, intent(out) ::  class_assgn(n)
   integer              ::  send_class_NN(k)
   real                 ::  p1, p2, posterior1, posterior2
   integer              ::  rank, p, chunk
   integer              ::  ierror,status(mpi_status_size)

   do test_row = 1,n

!Division on each block
   chunk=n/p
  ! if(rank < (n-chunk*p))then
  ! chunk=chunk+1
  ! endif
   print*,'reached knn'
    call distance_NN(train,test,n,m,test_row,k,send_class_NN,trainC,rank,p,chunk)

if(rank==0)then

   count_1 = 0
   count_2 = 0
   do i = 1, k
      if (send_class_NN(i)==1)then
        count_1 = count_1 + 1
      elseif (send_class_NN(i)==2)then
        count_2 = count_2 + 1
      endif
   enddo          
!ERROR CHECKING IF OCCURRED BY CHANCE    
   if ((count_1+count_2)==k) goto 1
   stop
 1 continue

!Posteriors
   posterior1 = FLOAT(count_1)*p1/(FLOAT(count_1)*p1 + FLOAT(count_2)*p2)
   posterior2 = FLOAT(count_2)*p2/(FLOAT(count_1)*p1 + FLOAT(count_2)*p2)
   
   if (posterior1>posterior2) then
      class_assgn(test_row) = 1
   endif
   if (posterior2>posterior1) then
      class_assgn(test_row) = 2
   endif
   if (posterior1==posterior2) then
      class_assgn(test_row) = 0
   endif
endif

  enddo

end subroutine nearest_neighbours
