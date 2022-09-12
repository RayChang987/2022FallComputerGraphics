#pragma once


class board {
private:
	static const int W = 8, H = 8;

	const int blue = 1, pink = 2, empty = 0;
public:
	void board_init();
	void check_board();
	void change_board(int x, int y, int val);
	int board[W][H];

};