      IntegerfunctionLENGTH(Needle)
      Character*(*)Needle
      LENGTH=LEN(Needle)
09999 IF(.NOT.(LENGTH.GT.0.AND.Needle(LENGTH:LENGTH).EQ.' '))GOTO09998
      LENGTH=LENGTH-1
      GOTO09999
09998 CONTINUE
      Return
      end
