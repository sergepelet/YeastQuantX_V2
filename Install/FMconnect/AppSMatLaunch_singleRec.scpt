FasdUAS 1.101.10   ��   ��    k             l     ����  O       	  r     
  
 b        b        m       �    [  n        4    �� 
�� 
ccel  m   	 
   �    E x p N u m  1    ��
�� 
pCRW  m       �    ]  o      ���� 
0 num Num 	 m       2                                                                                  FMP7  alis    �  Macintosh HD               ɟ��H+   Ve�FileMaker Pro Advanced.app                                      Vf����        ����  	                FileMaker Pro 11 Advanced     ɟ�v      ���x     Ve� J�  PMacintosh HD:Applications: FileMaker Pro 11 Advanced: FileMaker Pro Advanced.app  6  F i l e M a k e r   P r o   A d v a n c e d . a p p    M a c i n t o s h   H D  AApplications/FileMaker Pro 11 Advanced/FileMaker Pro Advanced.app   / ��  ��  ��        l     ��������  ��  ��        l    ����  r         m     ! ! � " " < N u m b e r   o f   r e c o r d   t o   a n a l y z e :   1   o      ���� 0 warnstr warnStr��  ��     # $ # l   " %���� % I   "�� & '
�� .sysodlogaskr        TEXT & o    ���� 0 warnstr warnStr ' �� ( )
�� 
btns ( J     * *  + , + m     - - � . .  S t o p ,  /�� / m     0 0 � 1 1  C o n t i n u e��   ) �� 2��
�� 
dflt 2 m    ���� ��  ��  ��   $  3 4 3 l  # 5 5���� 5 Z   # 5 6 7���� 6 =  # , 8 9 8 l  # ( :���� : n   # ( ; < ; 1   $ (��
�� 
bhit < l  # $ =���� = 1   # $��
�� 
rslt��  ��  ��  ��   9 m   ( + > > � ? ?  S t o p 7 L   / 1����  ��  ��  ��  ��   4  @ A @ l     ��������  ��  ��   A  B C B l     �� D E��   D Z Tdisplay dialog "Analysis Platform" buttons {"Mac", "PC", "cluster"} default button 1    E � F F � d i s p l a y   d i a l o g   " A n a l y s i s   P l a t f o r m "   b u t t o n s   { " M a c " ,   " P C " ,   " c l u s t e r " }   d e f a u l t   b u t t o n   1 C  G H G l     �� I J��   I ; 5set the Platform to the button returned of the result    J � K K j s e t   t h e   P l a t f o r m   t o   t h e   b u t t o n   r e t u r n e d   o f   t h e   r e s u l t H  L M L l     �� N O��   N O Idisplay dialog "Parallel Analysis" buttons {"Yes", "No"} default button 2    O � P P � d i s p l a y   d i a l o g   " P a r a l l e l   A n a l y s i s "   b u t t o n s   { " Y e s " ,   " N o " }   d e f a u l t   b u t t o n   2 M  Q R Q l     ��������  ��  ��   R  S T S l     �� U V��   U 8 2if the button returned of the result is "Yes" then    V � W W d i f   t h e   b u t t o n   r e t u r n e d   o f   t h e   r e s u l t   i s   " Y e s "   t h e n T  X Y X l     �� Z [��   Z ) #	set Platform to Platform & "_para"    [ � \ \ F 	 s e t   P l a t f o r m   t o   P l a t f o r m   &   " _ p a r a " Y  ] ^ ] l     �� _ `��   _  end if    ` � a a  e n d   i f ^  b c b l  6 Z d���� d O   6 Z e f e k   < Y g g  h i h r   < G j k j I  < C�� l��
�� .earsffdralis        afdr l m   < ?��
�� afdrcusr��   k o      ���� 0 
homefolder 
HomeFolder i  m n m r   H W o p o n   H S q r q 1   O S��
�� 
strq r l  H O s���� s n   H O t u t 1   K O��
�� 
psxp u o   H K���� 0 
homefolder 
HomeFolder��  ��   p o      ����  0 homefolderunix HomeFolderUnix n  v�� v l  X X�� w x��   w # display dialog HomeFolderUnix    x � y y : d i s p l a y   d i a l o g   H o m e F o l d e r U n i x��   f m   6 9 z z�                                                                                  MACS  alis    t  Macintosh HD               ɟ��H+   Jc
Finder.app                                                      J� �[��        ����  	                CoreServices    ɟ�v      �[ja     Jc JV JU  6Macintosh HD:System: Library: CoreServices: Finder.app   
 F i n d e r . a p p    M a c i n t o s h   H D  &System/Library/CoreServices/Finder.app  / ��  ��  ��   c  { | { l  [ h }���� } r   [ h ~  ~ b   [ d � � � b   [ b � � � o   [ ^����  0 homefolderunix HomeFolderUnix � m   ^ a � � � � � D L i b r a r y / F M c o n n e c t / Y Q l i n k . c o m m a n d     � o   b c���� 
0 num Num  o      ���� 0 	scriptstr 	scriptStr��  ��   |  � � � l     �� � ���   �  display dialog scriptStr    � � � � 0 d i s p l a y   d i a l o g   s c r i p t S t r �  � � � l  i w ����� � O   i w � � � I  o v�� ���
�� .coredoscnull��� ��� ctxt � o   o r���� 0 	scriptstr 	scriptStr��   � m   i l � ��                                                                                      @ alis    l  Macintosh HD               ɟ��H+   J�Terminal.app                                                    J>T� Y<        ����  	                	Utilities     ɟ�v      � =     J� J�  2Macintosh HD:Applications: Utilities: Terminal.app    T e r m i n a l . a p p    M a c i n t o s h   H D  #Applications/Utilities/Terminal.app   / ��  ��  ��   �  ��� � l     ��������  ��  ��  ��       �� � ���   � ��
�� .aevtoappnull  �   � **** � �� ����� � ���
�� .aevtoappnull  �   � **** � k     w � �   � �   � �  # � �  3 � �  b � �  { � �  �����  ��  ��   �   �   ����  �� !���� - 0���������� > z������������ ��� ���
�� 
pCRW
�� 
ccel�� 
0 num Num�� 0 warnstr warnStr
�� 
btns
�� 
dflt�� 
�� .sysodlogaskr        TEXT
�� 
rslt
�� 
bhit
�� afdrcusr
�� .earsffdralis        afdr�� 0 
homefolder 
HomeFolder
�� 
psxp
�� 
strq��  0 homefolderunix HomeFolderUnix�� 0 	scriptstr 	scriptStr
�� .coredoscnull��� ��� ctxt�� x� �*�,��/%�%E�UO�E�O����lv�l� O�a ,a   hY hOa  a j E` O_ a ,a ,E` OPUO_ a %�%E` Oa  	_ j U ascr  ��ޭ