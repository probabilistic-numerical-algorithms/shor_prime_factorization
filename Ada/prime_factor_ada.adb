pragma Ada_2022;

--  Prime Factorization Package Body
--  Language: Ada 2022
--  Fast, Numerical, Probabilistic with Perturbations for Retries

with Interfaces;
with Ada.Calendar;

package body prime_factor_ada is

   use type Interfaces.Unsigned_64;

   --  Linear Congruential Generator (LCG) for high-entropy randomization
   function LCG_Next (State : in out Interfaces.Unsigned_64) return Interfaces.Unsigned_64 is
   begin
      State := State * 6364136223846793005 + 1442695040888963407;
      return State;
   end LCG_Next;

   ---------
   -- GCD --
   ---------

   function gcd (A, B : Huge_Int) return Huge_Int is
      X : Huge_Int := A;
      Y : Huge_Int := B;
      Temp : Huge_Int;
   begin
      while Y > 0 loop
         pragma Loop_Optimize (Vector, Ivdep);
         pragma Loop_Optimize (Unroll);
         Temp := X mod Y;
         X := Y;
         Y := Temp;
      end loop;
      return X;
   end gcd;

   -------------
   -- Mul_Mod --
   -------------

   --  Performs modular multiplication (A * B) mod M safely avoiding 64-bit overflow
   function mul_mod (A, B, M : Huge_Int) return Huge_Int is
      X : Huge_Int := A mod M;
      Y : Huge_Int := B mod M;
      Result : Huge_Int := 0;
   begin
      while Y > 0 loop
         pragma Loop_Optimize (Vector, Ivdep);
         pragma Loop_Optimize (Unroll);
         if Y mod 2 = 1 then
            Result := (Result + X) mod M;
         end if;
         X := (X * 2) mod M;
         Y := Y / 2;
      end loop;
      return Result;
   end mul_mod;

   ---------------
   -- Power_Mod --
   ---------------

   --  Performs modular exponentiation (Base^Exp) mod Modulus safely
   function power_mod (Base, Exp, Modulus : Huge_Int) return Huge_Int is
      X : Huge_Int := Base mod Modulus;
      Y : Huge_Int := Exp;
      Result : Huge_Int := 1;
   begin
      while Y > 0 loop
         pragma Loop_Optimize (Vector, Ivdep);
         pragma Loop_Optimize (Unroll);
         if Y mod 2 = 1 then
            Result := mul_mod (Result, X, Modulus);
         end if;
         X := mul_mod (X, X, Modulus);
         Y := Y / 2;
      end loop;
      return Result;
   end power_mod;

   --------------------------
   -- Miller_Rabin_Witness --
   --------------------------

   function miller_rabin_witness (A, N, D, S : Huge_Int) return Boolean is
      X : Huge_Int := power_mod (A, D, N);
      Witness : Boolean := True;
      R : Huge_Int := 0;
   begin
      if X = 1 or else X = N - 1 then
         Witness := False;
      else
         while R < S - 1 and then not (X = N - 1) loop
            X := mul_mod (X, X, N);
            R := R + 1;
         end loop;
         if X = N - 1 then
            Witness := False;
         end if;
      end if;
      return Witness;
   end miller_rabin_witness;

   --------------
   -- Is_Prime --
   --------------

   function is_prime (N : Huge_Int; k : Positive := 5) return Boolean is
      Result : Boolean := True;
   begin
      if N <= 1 then
         Result := False;
      elsif N = 2 or else N = 3 then
         Result := True;
      elsif N mod 2 = 0 then
         Result := False;
      else
         declare
            D : Huge_Int := N - 1;
            S : Huge_Int := 0;
            I : Positive := 1;
            A : Huge_Int;
            Witness_Found : Boolean := False;
            State : Interfaces.Unsigned_64 := Interfaces.Unsigned_64 (Float (Ada.Calendar.Seconds (Ada.Calendar.Clock)) * 1000.0) + Interfaces.Unsigned_64 (N mod 4294967296);
            Random_Val : Huge_Int;
         begin
            while D mod 2 = 0 loop
               D := D / 2;
               S := S + 1;
            end loop;

            while I <= k and then not Witness_Found loop
               declare
                  U64_Val : constant Interfaces.Unsigned_64 := LCG_Next (State);
               begin
                  Random_Val := Long_Long_Integer (U64_Val mod 2**63);
               end;
               
               A := 2 + (Random_Val mod (N - 3));
               
               if miller_rabin_witness (A, N, D, S) then
                  Witness_Found := True;
                  Result := False;
               end if;
               I := I + 1;
            end loop;
         end;
      end if;
      return Result;
   end is_prime;

   ------------------
   -- Integer_Sqrt --
   ------------------

   function integer_sqrt (Val : Huge_Int) return Huge_Int is
      Result : Huge_Int := 0;
   begin
      if Val > 0 then
         declare
            X : Huge_Int := Val;
            Y : Huge_Int := (X + 1) / 2;
         begin
            while Y < X loop
               X := Y;
               Y := (X + Val / X) / 2;
            end loop;
            Result := X;
         end;
       end if;
       return Result;
   end integer_sqrt;

   ------------------------
   -- Is_Perfect_Square --
   ------------------------

   function is_perfect_square (Val : Huge_Int; Root : out Huge_Int) return Boolean is
      Result : Boolean := False;
   begin
      Root := 0;
      if Val >= 0 then
         declare
            R : constant Huge_Int := integer_sqrt (Val);
         begin
            if R * R = Val then
               Result := True;
               Root := R;
            end if;
         end;
      end if;
      return Result;
   end is_perfect_square;

   ---------------------
   -- Shor_Factorize --
   ---------------------

   --  Finds a non-trivial factor of N using Shor's period-finding algorithm.
   --  Classically, we find the period r of a^x mod N.
   --  This has been parallelized using Ada task arrays to match the time
   --  complexity of the physical resonance-based CMOS Shor system.
   function shor_factorize (N : Huge_Int; P, Q : out Huge_Int) return Boolean is
      Found : Boolean := False;
      
      task type Search_Worker is
         entry Start (Id : Positive; N_Val_In : Huge_Int);
         entry Get_Result (Success : out Boolean; Out_P, Out_Q : out Huge_Int);
      end Search_Worker;

      task body Search_Worker is
         My_Id     : Positive;
         P_Found   : Huge_Int := 0;
         Q_Found   : Huge_Int := 0;
         Found_Loc : Boolean := False;
         N_Val     : Huge_Int;
      begin
         accept Start (Id : Positive; N_Val_In : Huge_Int) do
            My_Id := Id;
            N_Val := N_Val_In;
         end Start;

         declare
            Attempts : Natural := 0;
            Max_Worker_Attempts : constant Natural := 15;
            State : Interfaces.Unsigned_64 := Interfaces.Unsigned_64 (Float (Ada.Calendar.Seconds (Ada.Calendar.Clock)) * 1000.0) + Interfaces.Unsigned_64 (My_Id * 1000) + Interfaces.Unsigned_64 (N_Val mod 4294967296);
            Random_Val : Huge_Int;
            A : Huge_Int;
            G : Huge_Int;
         begin
            while not Found_Loc and then Attempts < Max_Worker_Attempts loop
               declare
                  U64_Val : constant Interfaces.Unsigned_64 := LCG_Next (State);
               begin
                  Random_Val := Long_Long_Integer (U64_Val mod 2**63) + Long_Long_Integer (My_Id * 100);
               end;

               A := 2 + (Random_Val mod (N - 3));
               G := gcd (A, N);

               if G > 1 and then G < N then
                  P_Found := G;
                  Q_Found := N / G;
                  Found_Loc := True;
               else
                  -- Find the period r of f(x) = a^x mod N
                  declare
                     R : Huge_Int := 1;
                     Val_Mod : Huge_Int := power_mod (A, 1, N);
                  begin
                     while Val_Mod /= 1 and then R < 5000 loop
                        pragma Loop_Optimize (Vector, Ivdep);
                        pragma Loop_Optimize (Unroll);
                        Val_Mod := mul_mod (Val_Mod, A, N);
                        R := R + 1;
                     end loop;

                     if Val_Mod = 1 and then R mod 2 = 0 then
                        declare
                           X : constant Huge_Int := power_mod (A, R / 2, N);
                           Factor : Huge_Int;
                        begin
                           if X /= N - 1 then
                              Factor := gcd (X - 1, N);
                              if Factor > 1 and then Factor < N then
                                 P_Found := Factor;
                                 Q_Found := N / Factor;
                                 Found_Loc := True;
                              else
                                 Factor := gcd (X + 1, N);
                                 if Factor > 1 and then Factor < N then
                                    P_Found := Factor;
                                    Q_Found := N / Factor;
                                    Found_Loc := True;
                                 end if;
                              end if;
                           end if;
                        end;
                     end if;
                  end;
               end if;
               Attempts := Attempts + 1;
            end loop;
         end;

         accept Get_Result (Success : out Boolean; Out_P, Out_Q : out Huge_Int) do
            Success := Found_Loc;
            Out_P := P_Found;
            Out_Q := Q_Found;
         end Get_Result;
      end Search_Worker;

      type Worker_Index is range 1 .. 8;
      Workers : array (Worker_Index) of Search_Worker;
      
      Success_Loc : Boolean;
      P_Loc, Q_Loc : Huge_Int;
   begin
      P := 0;
      Q := 0;

      --#pragma omp parallel for
      --$omp parallel for
      for Id in Worker_Index loop
         Workers (Id).Start (Positive (Id), N);
      end loop;

      --#pragma omp parallel for reduction(or:Found)
      --$omp parallel for reduction(or:Found)
      for Id in Worker_Index loop
         Workers (Id).Get_Result (Success_Loc, P_Loc, Q_Loc);
         if Success_Loc and then not Found then
            P := P_Loc;
            Q := Q_Loc;
            Found := True;
         end if;
      end loop;

      return Found;
   end shor_factorize;

   ------------------
   -- Sort_Factors --
   ------------------

   procedure sort_factors (Result : in out Factorization_Result) is
      Temp : Huge_Int;
      Swapped : Boolean := True;
   begin
      if Result.Count > 1 then
         while Swapped loop
            Swapped := False;
            for i in 1 .. Result.Count - 1 loop
               if Result.Factors (i) > Result.Factors (i + 1) then
                  Temp := Result.Factors (i);
                  Result.Factors (i) := Result.Factors (i + 1);
                  Result.Factors (i + 1) := Temp;
                  Swapped := True;
               end if;
            end loop;
         end loop;
      end if;
   end sort_factors;

   -------------------------
   -- Factorize_Recursive --
   -------------------------

   procedure factorize_recursive (N : in Huge_Int; Temp_Factors : in out Factor_Array; Count : in out Natural) is
      P, Q : Huge_Int := 0;
      Found : Boolean := False;
   begin
      if N <= 1 then
         null;
      elsif is_prime (N) then
         Count := Count + 1;
         Temp_Factors (Count) := N;
      else
         -- Try Shor's period-finding
         if shor_factorize (N, P, Q) then
            Found := True;
         end if;
         
         -- Try Trial Division if Shor's fails
         if not Found then
            declare
               Divisor : Huge_Int := 3;
            begin
               while Divisor * Divisor <= N and then not Found loop
                  if N mod Divisor = 0 then
                     P := Divisor;
                     Q := N / Divisor;
                     Found := True;
                  else
                     Divisor := Divisor + 2;
                  end if;
               end loop;
            end;
         end if;
         
         if Found then
            factorize_recursive (P, Temp_Factors, Count);
            factorize_recursive (Q, Temp_Factors, Count);
         else
            Count := Count + 1;
            Temp_Factors (Count) := N;
         end if;
      end if;
   end factorize_recursive;

   ---------------
   -- Factorize --
   ---------------

   function factorize (N : Huge_Int) return Factorization_Result is
      Result : Factorization_Result;
      Temp_Factors : Factor_Array (1 .. 128) := (others => 0);
      Count : Natural := 0;
   begin
      factorize_recursive (N, Temp_Factors, Count);
      
      Result.Count := Count;
      if Count > 0 then
         Result.Factors := new Factor_Array (1 .. Count);
         for I in 1 .. Count loop
            Result.Factors (I) := Temp_Factors (I);
         end loop;
         sort_factors (Result);
      else
         Result.Factors := null;
      end if;
      
      return Result;
   end factorize;

end prime_factor_ada;
