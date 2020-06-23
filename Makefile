all: paddle t1 t2 t3 t4 t5 t6 t7 t8 t9

paddle: paddle.cpp
	g++ paddle.cpp -Wall -opaddle -lX11

t1: t1.cpp
	g++ t1.cpp -Wall -ot1 -lX11

t2: t2.cpp
	g++ t2.cpp -Wall -ot2 -lX11

t3: t3.cpp
	g++ t3.cpp -Wall -ot3 -lX11

t4: t4.cpp
	g++ t4.cpp -Wall -ot4 -lX11

t5: t5.cpp
	g++ t5.cpp -Wall -ot5 -lX11

t6: t6.cpp
	g++ t6.cpp -Wall -ot6 -lX11

t7: t7.cpp
	g++ t7.cpp -Wall -ot7 -lX11

t8: t8.cpp
	g++ t8.cpp -Wall -ot8 -lX11

t9: t9.cpp
	g++ t9.cpp -Wall -ot9 -lX11

clean:
	rm -f paddle
	rm -f t1
	rm -f t2
	rm -f t3
	rm -f t4
	rm -f t5
	rm -f t6
	rm -f t7
	rm -f t8
	rm -f t9
