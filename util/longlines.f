	character*80 line

10	read(5,'(A)',END=20) line
	write(6,'(A)') line
	go to 10
20	stop
	end
