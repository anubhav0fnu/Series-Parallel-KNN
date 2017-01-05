subroutine nearest_neighbours(train,test,n,m,k,class_assgn,trainC,p1,p2)
   implicit none
   integer              ::  test_row, i, k, n, m, count_1, count_2
   real   , intent(in)  ::  train(n,m), test(n,m)
   integer, intent(in)  ::  trainC(n)
   integer, intent(out) ::  class_assgn(n)
   integer              ::  class_NN(k)
   real                 ::  p1, p2, posterior1, posterior2

   do test_row = 1,n
     call distance_NN(train,test,n,m,test_row,k,class_NN,trainC)

   count_1 = 0
   count_2 = 0
   do i = 1, k
      if (class_NN(i)==1)then
        count_1 = count_1 + 1
      elseif (class_NN(i)==2)then
        count_2 = count_2 + 1
      endif
   enddo          
    
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
  enddo

end subroutine nearest_neighbours
