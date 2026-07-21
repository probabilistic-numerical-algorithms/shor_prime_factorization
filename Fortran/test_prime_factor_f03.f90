!  Prime Factorization Test Suite
!  Language: Fortran 2003

program test_prime_factor_f03
  use prime_factor_f03
  implicit none

  integer, parameter :: test_cases_count = 11
  integer(kind=8) :: test_numbers(test_cases_count)
  integer(kind=8), allocatable :: factors(:)
  integer :: count, i, j, pass_count, test_count
  integer(kind=8) :: product
  logical :: all_prime

  ! Define a range of test numbers to factorize
  test_numbers = (/ 2_8, 15_8, 28_8, 101_8, 121_8, 1043_8, 840_8, 10007_8, 65535_8, 1000000_8, 4294967297_8 /)

  pass_count = 0
  test_count = 0

  write(*,*) "============================================================="
  write(*,*) "  PRIME FACTORIZATION COMPREHENSIVE TEST SUITE (F2003)"
  write(*,*) "============================================================="
  write(*,*) ""

  do i = 1, test_cases_count
     test_count = test_count + 1
     write(*,*) "-------------------------------------------------------------"
     write(*, "(A,I0,A,I0)") "TEST ", i, ": Factorizing N = ", test_numbers(i)
     write(*,*) "-------------------------------------------------------------"

     call factorize(test_numbers(i), factors)
     count = size(factors)

     ! Verify factors
     product = 1_8
     all_prime = .true.

     write(*, fmt='(A)', advance='no') "Factors found: "
     if (count == 0) then
        write(*, fmt='(A)') "(none)"
     else
        do j = 1, count
           write(*, fmt='(I0,A)', advance='no') factors(j), " "
           product = product * factors(j)
           if (.not. is_prime(factors(j))) then
              all_prime = .false.
           end if
        end do
        write(*,*) ""
     end if

     write(*, "(A,I0)") "Product of factors: ", product
     write(*, "(A,L1)") "Are all factors prime? ", all_prime

     if (product == test_numbers(i) .and. all_prime) then
        write(*,*) "STATUS: PASS"
        pass_count = pass_count + 1
     else
        write(*,*) "STATUS: FAIL"
        if (product /= test_numbers(i)) then
           write(*,*) "  Reason: Product of factors does not equal original number!"
        end if
        if (.not. all_prime) then
           write(*,*) "  Reason: One or more factors are composite!"
        end if
      end if
      write(*,*) ""
  end do

  write(*,*) "============================================================="
  write(*,*) "  TEST SUITE SUMMARY"
  write(*,*) "============================================================="
  write(*,*) "Tests run:   ", test_count
  write(*,*) "Tests passed:", pass_count

  if (pass_count == test_count) then
     write(*,*) "OVERALL STATUS: ALL TESTS PASSED"
  else
     write(*,*) "OVERALL STATUS: FAILURE"
  end if
  write(*,*) "============================================================="

end program test_prime_factor_f03
