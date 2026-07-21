pragma Ada_2022;

with ada.text_io; use ada.text_io;
with prime_factor_ada; use prime_factor_ada;

procedure test_prime_factor_ada is

   type Test_Case is record
      N : Huge_Int;
      Expected_Count : Natural;
   end record;

   --  A selection of test cases
   Test_Suite : constant array (1 .. 16) of Test_Case := (
      (N => 4, Expected_Count => 2),           -- 2 * 2
      (N => 12, Expected_Count => 3),          -- 2 * 2 * 3
      (N => 15, Expected_Count => 2),          -- 3 * 5
      (N => 35, Expected_Count => 2),          -- 5 * 7
      (N => 143, Expected_Count => 2),         -- 11 * 13
      (N => 1024, Expected_Count => 10),       -- 2^10
      (N => 4087, Expected_Count => 2),        -- 61 * 67
      (N => 999999, Expected_Count => 7),      -- 3 * 3 * 3 * 7 * 11 * 13 * 37
      (N => 123456789, Expected_Count => 4),   -- 3 * 3 * 3607 * 3803
      (N => 11111111111, Expected_Count => 2), -- 21649 * 513239
      (N => 17, Expected_Count => 1),          -- prime 17
      (N => 1_152_921_504_606_641_767,
       Expected_Count => 2),                   -- 1_073_741_371 * 1_073_742_277
      (N => 2_305_843_009_212_522_371,
       Expected_Count => 2),                   -- 1_518_494_123 * 1_518_494_123
      (N => 4_611_686_018_425_746_943,
       Expected_Count => 2),                   -- 2_147_482_367 * 2_147_484_929
      (N => 4_611_686_018_425_746_947,
       Expected_Count => 3),                   -- 109 *  181 *  233751635583443
      (N => 4_611_686_018_425_746_979,
       Expected_Count => 5)                    -- 19 *  23 *  397 *  910051 *  29209361
   );

   Result : Factorization_Result;
   Product : Huge_Int;
   All_Prime : Boolean;
   Pass_Count : Natural := 0;
   Test_Count : Natural := 0;

begin
   put_line ("=============================================================");
   put_line ("  PRIME FACTORIZATION COMPREHENSIVE TEST SUITE (ADA 2022)");
   put_line ("=============================================================");
   new_line;

   for i in Test_Suite'range loop
      Test_Count := Test_Count + 1;
      put_line ("-------------------------------------------------------------");
      put_line ("TEST " & i'image & ": Factorizing N = " & Test_Suite (i).N'Image);
      put_line ("-------------------------------------------------------------");

      Result := factorize (Test_Suite (i).N);
      
      -- Verify the result
      Product := 1;
      All_Prime := True;
      
      put ("Factors found: ");
      if Result.Count = 0 then
         put ("(none)");
      else
         for j in 1 .. Result.Count loop
            put (" " & Result.Factors (j)'Image);
            if j < Result.Count then
               put (" *");
            end if;
            Product := Product * Result.Factors (j);
            if not is_prime (Result.Factors (j)) then
               All_Prime := False;
            end if;
         end loop;
      end if;
      new_line;

      put_line ("Product of factors: " & Product'Image);
      put_line ("Are all factors prime? " & All_Prime'image);

      if Product = Test_Suite (i).N and then All_Prime then
         put_line ("STATUS: PASS");
         Pass_Count := Pass_Count + 1;
      else
         put_line ("STATUS: FAIL");
         if Product /= Test_Suite (i).N then
            put_line ("  Reason: Product of factors does not equal original number!");
         end if;
         if not All_Prime then
            put_line ("  Reason: One or more factors are composite!");
         end if;
      end if;
      new_line;
   end loop;

   put_line ("=============================================================");
   put_line ("  TEST SUITE SUMMARY");
   put_line ("=============================================================");
   put_line ("Tests run:   " & Test_Count'image);
   put_line ("Tests passed: " & Pass_Count'image);
   
   if Pass_Count = Test_Count then
      put_line ("OVERALL STATUS: ALL TESTS PASSED");
   else
      put_line ("OVERALL STATUS: FAILURE");
   end if;
   put_line ("=============================================================");

end test_prime_factor_ada;
