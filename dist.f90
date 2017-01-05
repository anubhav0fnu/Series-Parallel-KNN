!FIRST Approach
subroutine distance_NN(train,test,n,m,test_row,k,class_NN,trainC)
   implicit none
   integer                :: i, j, k, n, m, train_row, test_row, col
   real,intent(IN)        :: train(n,m), test(n,m)
   integer,intent(IN)     :: trainC(n)
   integer,intent(out)    :: class_NN(k)
   real                   :: dist(k),d 
!initialise the arrays
   do i  = 1, k
      dist(i)     = 1.E+15
      class_NN(i) = 0
   enddo
   do train_row = 1,n
      d = 0.0
      do col = 1, m
         d = d + (train(train_row,col) - test(test_row,col))**2
      enddo
      do i = 1, k
         if (d <  dist(i)) then
            do j = k,i+1,-1
               if (j <= k) then
                   dist(j)      = dist(j-1)
                   class_NN(j)  = class_NN(j-1)
               endif
            enddo
         dist(i)      = d
         class_NN(i)  = trainC(train_row)
         goto 1
         endif
      enddo    
 1       continue
   enddo

end subroutine distance_NN
