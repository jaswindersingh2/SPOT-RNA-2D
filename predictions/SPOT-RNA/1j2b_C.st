#Name: 1j2b_C
#Length:  77 
#PageNumber: 2
GGGCCCGUGGUCUAGUUGGUCAUGACGCCGCCCUUACGAGGCGGAGGUCCGGGGUUCAAGUCCCCGCGGGCCCACCA
((((((([.((((][.......))))((((((.......))))))...]((((((...)..))))))))))))....
SSSSSSSMMSSSSHHHHHHHHHSSSSSSSSSSHHHHHHHSSSSSSMMMMSSSSSSHHHSBBSSSSSSSSSSSSEEEE
NNNNNNNKNNNNNKKNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNKNNNNNNNNNNNNNNNNNNNNNNNNNNNN
S1 1..7 "GGGCCCG" 67..73 "CGGGCCC"
S2 10..13 "GUCU" 23..26 "UGAC"
S3 27..32 "GCCGCC" 40..45 "GGCGGA"
S4 50..54 "CGGGG" 62..66 "CCCCG"
S5 55..55 "U" 59..59 "A"
H1 14..22 "AGUUGGUCA" (13,23) U:U PK{1,2}
H2 33..39 "CUUACGA" (32,40) C:G 
H3 56..58 "UCA" (55,59) U:A 
B1 60..61 "GU" (59,55) A:U (62,54) C:G 
M1.1 8..9 "UG" (7,67) G:C (10,26) G:C PK{1}
M1.2 27..26 "" (26,10) C:G (27,45) G:A 
M1.3 46..49 "GGUC" (45,27) A:G (50,66) C:G PK{2}
M1.4 67..66 "" (66,50) G:C (67,7) C:G 
E1 74..77 "ACCA" 
PK1 1bp 8..8 14..14 M1.1 8..9 H1 14..22
PK2 1bp 15..15 49..49 H1 14..22 M1.3 46..49
PK1.1 8 U 14 A
PK2.1 15 G 49 C
NCBP1 13 U 23 U S2
NCBP2 27 G 45 A S3
segment1 7bp 1..7 GGGCCCG 67..73 CGGGCCC
segment2 4bp 10..13 GUCU 23..26 UGAC
segment3 6bp 27..32 GCCGCC 40..45 GGCGGA
segment4 6bp 50..55 CGGGGU 59..66 AGUCCCCG