#include "board.h";
#include <iostream>
void board::board_init() {
	for (int i = 0; i < W; ++i) {
		for (int j = 0; j < H; ++j) {
			board[i][j] = empty;
		}
	}
	board[3][3] = board[4][4] = pink;
	board[3][4] = board[4][3] = blue;
}
void board::check_board() {
	for (int i = 0; i < W; ++i) {
		for (int j = 0; j < H; ++j) {
			std::cout << board[i][j] << " ";
		}
		std::cout << std::endl;
	}
}
void board::change_board(int x, int y, int val) {
	board[x][y] = val;
}