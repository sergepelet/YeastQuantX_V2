FasdUAS 1.101.10   ��   ��    k             l     ��������  ��  ��        l    D 	���� 	 O     D 
  
 k    C       r        l   	 ����  I   	�� ��
�� .corecnte****       ****  m    ��
�� 
crow��  ��  ��    o      ���� 0 foundrec FoundRec   ��  Y    C ��  ��  k    >       I   �� ��
�� .FMPRGOTOnull���    obj   4    �� 
�� 
crow  o    ���� 0 i  ��     ��  Z    >  ��   =   "   !   o     ���� 0 i   ! m     !����   r   % / " # " b   % - $ % $ m   % & & & � ' '  [ % n   & , ( ) ( 4   ) ,�� *
�� 
ccel * m   * + + + � , ,  E x p N u m ) 1   & )��
�� 
pCRW # o      ���� 
0 num Num��    r   2 > - . - b   2 < / 0 / b   2 5 1 2 1 o   2 3���� 
0 num Num 2 m   3 4 3 3 � 4 4  , 0 n   5 ; 5 6 5 4   8 ;�� 7
�� 
ccel 7 m   9 : 8 8 � 9 9  E x p N u m 6 1   5 8��
�� 
pCRW . o      ���� 
0 num Num��  �� 0 i    m    ����   o    ���� 0 foundrec FoundRec��  ��    m      : :2                                                                                  FMP7  alis    �  Macintosh HD               ɟ��H+   Ve�FileMaker Pro Advanced.app                                      Vf����        ����  	                FileMaker Pro 11 Advanced     ɟ�v      ���x     Ve� J�  PMacintosh HD:Applications: FileMaker Pro 11 Advanced: FileMaker Pro Advanced.app  6  F i l e M a k e r   P r o   A d v a n c e d . a p p    M a c i n t o s h   H D  AApplications/FileMaker Pro 11 Advanced/FileMaker Pro Advanced.app   / ��  ��  ��     ; < ; l  E J =���� = r   E J > ? > b   E H @ A @ o   E F���� 
0 num Num A m   F G B B � C C  ] ? o      ���� 
0 num Num��  ��   <  D E D l     ��������  ��  ��   E  F G F l  K P H���� H r   K P I J I b   K N K L K m   K L M M � N N : N u m b e r   o f   r e c o r d   t o   a n a l y z e :   L o   L M���� 0 foundrec FoundRec J o      ���� 0 warnstr warnStr��  ��   G  O P O l  Q e Q���� Q I  Q e�� R S
�� .sysodlogaskr        TEXT R o   Q R���� 0 warnstr warnStr S �� T U
�� 
btns T J   S [ V V  W X W m   S V Y Y � Z Z  S t o p X  [�� [ m   V Y \ \ � ] ]  C o n t i n u e��   U �� ^��
�� 
dflt ^ m   ^ _���� ��  ��  ��   P  _ ` _ l  f z a���� a Z   f z b c���� b =  f q d e d l  f m f���� f n   f m g h g 1   i m��
�� 
bhit h l  f i i���� i 1   f i��
�� 
rslt��  ��  ��  ��   e m   m p j j � k k  S t o p c L   t v����  ��  ��  ��  ��   `  l m l l     ��������  ��  ��   m  n o n l     �� p q��   p Z Tdisplay dialog "Analysis Platform" buttons {"Mac", "PC", "cluster"} default button 1    q � r r � d i s p l a y   d i a l o g   " A n a l y s i s   P l a t f o r m "   b u t t o n s   { " M a c " ,   " P C " ,   " c l u s t e r " }   d e f a u l t   b u t t o n   1 o  s t s l     �� u v��   u ; 5set the Platform to the button returned of the result    v � w w j s e t   t h e   P l a t f o r m   t o   t h e   b u t t o n   r e t u r n e d   o f   t h e   r e s u l t t  x y x l     �� z {��   z O Idisplay dialog "Parallel Analysis" buttons {"Yes", "No"} default button 2    { � | | � d i s p l a y   d i a l o g   " P a r a l l e l   A n a l y s i s "   b u t t o n s   { " Y e s " ,   " N o " }   d e f a u l t   b u t t o n   2 y  } ~ } l     ��������  ��  ��   ~   �  l     �� � ���   � 8 2if the button returned of the result is "Yes" then    � � � � d i f   t h e   b u t t o n   r e t u r n e d   o f   t h e   r e s u l t   i s   " Y e s "   t h e n �  � � � l     �� � ���   � ) #	set Platform to Platform & "_para"    � � � � F 	 s e t   P l a t f o r m   t o   P l a t f o r m   &   " _ p a r a " �  � � � l     �� � ���   �  end if    � � � �  e n d   i f �  � � � l  { � ����� � O   { � � � � k   � � � �  � � � r   � � � � � I  � ��� ���
�� .earsffdralis        afdr � m   � ���
�� afdrcusr��   � o      ���� 0 
homefolder 
HomeFolder �  � � � r   � � � � � n   � � � � � 1   � ���
�� 
strq � l  � � ����� � n   � � � � � 1   � ���
�� 
psxp � o   � ����� 0 
homefolder 
HomeFolder��  ��   � o      ����  0 homefolderunix HomeFolderUnix �  ��� � l  � ��� � ���   � # display dialog HomeFolderUnix    � � � � : d i s p l a y   d i a l o g   H o m e F o l d e r U n i x��   � m   { ~ � ��                                                                                  MACS  alis    t  Macintosh HD               ɟ��H+   Jc
Finder.app                                                      J� �[��        ����  	                CoreServices    ɟ�v      �[ja     Jc JV JU  6Macintosh HD:System: Library: CoreServices: Finder.app   
 F i n d e r . a p p    M a c i n t o s h   H D  &System/Library/CoreServices/Finder.app  / ��  ��  ��   �  � � � l  � � ����� � r   � � � � � b   � � � � � b   � � � � � o   � �����  0 homefolderunix HomeFolderUnix � m   � � � � � � � D L i b r a r y / F M c o n n e c t / Y Q l i n k . c o m m a n d     � o   � ����� 
0 num Num � o      ���� 0 	scriptstr 	scriptStr��  ��   �  � � � l     �� � ���   �  display dialog scriptStr    � � � � 0 d i s p l a y   d i a l o g   s c r i p t S t r �  � � � l  � � ����� � O   � � � � � I  � ��� ���
�� .coredoscnull��� ��� ctxt � o   � ����� 0 	scriptstr 	scriptStr��   � m   � � � ��                                                                                      @ alis    l  Macintosh HD               ɟ��H+   J�Terminal.app                                                    J>T� Y<        ����  	                	Utilities     ɟ�v      � =     J� J�  2Macintosh HD:Applications: Utilities: Terminal.app    T e r m i n a l . a p p    M a c i n t o s h   H D  #Applications/Utilities/Terminal.app   / ��  ��  ��   �  � � � l     ��������  ��  ��   �  ��� � l     ��������  ��  ��  ��       �� � ���   � ��
�� .aevtoappnull  �   � **** � �� ����� � ���
�� .aevtoappnull  �   � **** � k     � � �   � �  ; � �  F � �  O � �  _ � �  � � �  � � �  �����  ��  ��   � ���� 0 i   � # :�������� &���� +�� 3 8 B M���� Y \��������� j ��~�}�|�{�z�y ��x ��w
�� 
crow
�� .corecnte****       ****�� 0 foundrec FoundRec
�� .FMPRGOTOnull���    obj 
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
� 
bhit
�~ afdrcusr
�} .earsffdralis        afdr�| 0 
homefolder 
HomeFolder
�{ 
psxp
�z 
strq�y  0 homefolderunix HomeFolderUnix�x 0 	scriptstr 	scriptStr
�w .coredoscnull��� ��� ctxt�� �� A�j E�O 6k�kh  *�/j O�k  �*�,��/%E�Y ��%*�,��/%E�[OY��UO��%E�O��%E�O��a a lva la  O_ a ,a   hY hOa  a j E` O_ a ,a ,E` OPUO_ a %�%E`  Oa ! 	_  j "Uascr  ��ޭ