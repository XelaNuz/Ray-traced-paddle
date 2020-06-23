all: paddle

paddle: paddle.cpp
	g++ paddle.cpp -Wall -opaddle -lX11

clean:
	rm -f paddle
