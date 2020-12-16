      SUBROUTINEVTIME(XTIME)
      REALXTIME
      externalsecond
      XTIME=SECOND(0.0)
      RETURN
      END
      realfunctionsecond(t)
      real(4)t,tcpu
      callcpu_time(tcpu)
      second=tcpu
      end
      SubroutineTremain(RTime,BegTime)
      Parameter(RunTime=86400.)
      RealXtime,Rtime,BegTime
      callVtime(Xtime)
      Rtime=RunTime-(Xtime-BegTime)
      Return
      End
      SubroutineDateTime(Date,Time)
      Character*(*)Date,Time
      IntegerHour,minute,second,hundredth
      Integer*4Month,Day,Year
      INTEGERDATE_TIME(8)
      CHARACTER(LEN=12)REAL_CLOCK(3)
      CHARACTER(10)t
      CHARACTER(5)z
      CALLDATE_AND_TIME(REAL_CLOCK(1),REAL_CLOCK(2),REAL_CLOCK(3),DATE_TIME)
      write(*,*)' REAL_CLOCK(1)=',REAL_CLOCK(1)
      write(*,*)' REAL_CLOCK(2)=',REAL_CLOCK(2)
      write(*,*)' REAL_CLOCK(3)=',REAL_CLOCK(3)
      write(*,*)' DATE_TIME=',DATE_TIME
      CALLDATE_AND_TIME(TIME=t,ZONE=z)
      write(*,*)' time=',t
      write(*,*)' zone=',z
      write(Date,'(I4,''/'',I2,''/'',I2)')Date_time(1),Date_time(2),Date_time(3)
      Hour=Date_time(5)
      minute=Date_time(6)
      second=Date_time(7)
      write(Time,'(I2,'':'',I2,'':'',I2)')Hour,minute,second
      Return
      End
      SubroutineBasTat(string)
      Character*(*)string
      Return
      End
      SubroutineDatTim(d,t)
      character*(*)d,t
      Return
      End
      FunctionLtime(itim)
      Ltime=0
      Return
      End
