#pragma once
#include "Loader.h"
#include<iostream>
#include <vector>
sf::Vector2i translator(sf::Vector2i point, const int dem) {
	return sf::Vector2i(point.x / dem, point.y / dem);
}
 void game::read_board_and_react() {
	for (int i = 0; i < gridN; ++i) {
		for (int j = 0; j < gridN; ++j) {
			if (main_board.board[i][j] == 0) { continue; }
			if (main_board.board[i][j] == 1) {
				window.draw(blue_figs[i * 8 + j]);
			}
			if (main_board.board[i][j] == 2) {
				window.draw(pink_figs[i * 8 + j]);
			}
		}
	}
}
const int dirX[8] = { 0, 0, 1, -1, 1, -1, 1, -1 };
const int dirY[8] = { 1, -1, 0, 0, 1, -1, -1, 1 };
void game::click_handler(const sf::Vector2i& point) {
	sf::Vector2i p = translator(point, winW/gridN);
	if (p.x >= 8 || p.x < 0 || p.y>=8 || p.y < 0) { std::cerr << "ERROR: Click Out OF Range!" << std::endl; return; }
	if (sta.turn ==sta.blue) {
		if (main_board.board[p.x][p.y] == 0) {
			bool legal = false;
			for (int i = 0; i < 8; ++i) {
				bool fnd_pink = false;
				bool flip = false;
				for (int j = 1; j < 8; ++j) {
					int nX = p.x + dirX[i] * j, nY = p.y + dirY[i] * j;
					if (nX < 0 || nX >= 8 || nY < 0 || nY >= 8) { break; }
					if (main_board.board[nX][nY] == sta.empty || (fnd_pink==0&&main_board.board[nX][nY] == sta.blue)) { break; }
					if(main_board.board[nX][nY]==sta.blue&&fnd_pink){
						legal = true;
						flip = true;
					}
					if (main_board.board[nX][nY] == sta.pink) { fnd_pink = true;}
				}
				if (flip) {
					for (int j = 1; j < 8; ++j) {
						int nX = p.x + dirX[i] * j, nY = p.y + dirY[i] * j;
						if (main_board.board[nX][nY] == sta.blue||main_board.board[nX][nY]==sta.empty) { break; }
						if (main_board.board[nX][nY] == sta.pink) { main_board.board[nX][nY] = sta.blue; sta.blue_q++; sta.pink_q--; }
					}
				}
			}
			if (legal) {
				main_board.board[p.x][p.y] = sta.blue;
				sound.play();
				sta.blue_q++;
				sta.change();
			}
			else { std::cout << "It's not a legal move!" << std::endl; }
		}
		
	}
	if (sta.turn == sta.pink) {
		if (main_board.board[p.x][p.y] == 0) {
			bool legal = false;
			for (int i = 0; i < 8; ++i) {
				bool fnd_blue = false;
				bool flip = false;
				for (int j = 1; j < 8; ++j) {
					int nX = p.x + dirX[i] * j, nY = p.y + dirY[i] * j;
					if (nX < 0 || nX >= 8 || nY < 0 || nY >= 8) { break; }
					if (main_board.board[nX][nY] == sta.empty || (fnd_blue == 0 && main_board.board[nX][nY] == sta.pink)) { break; }
					if (main_board.board[nX][nY] == sta.pink && fnd_blue) {
						legal = true;
						flip = true;
					}
					if (main_board.board[nX][nY] == sta.pink || main_board.board[nX][nY] == sta.empty) { break; }
					if (main_board.board[nX][nY] == sta.blue) { fnd_blue = true; }
				}
				if (flip) {
					for (int j = 1; j < 8; ++j) {
						int nX = p.x + dirX[i] * j, nY = p.y + dirY[i] * j;
						if (main_board.board[nX][nY] == sta.pink || main_board.board[nX][nY] == sta.empty) { break; }
						if (main_board.board[nX][nY] == sta.blue) { main_board.board[nX][nY] = sta.pink; sta.pink_q++; sta.blue_q--; }
					}
				}
			}
			if (legal) {
				main_board.board[p.x][p.y] = sta.pink;
				sta.pink_q++;
				sta.change();
				sound.play();
			}
			else { std::cout << "It's not a legal move!"<<std::endl; }
		}

	}
	sta.end_check();
	//std::cout << sta.blue_q << " " << sta.pink_q << std:: endl;
}
void game::refresh() {
	window.draw(background);
	read_board_and_react();
	window.display();
}
void game::fill_position_map() {
	int x = 4, y = 6; //35 for winH = winW = 560, gridN = 80
	for (int i = 0; i < gridN; ++i) {
		for (int j = 0; j < gridN; ++j) {
			position_map[i][j].x = x + i * 70;
			position_map[i][j].y = y + j * 70;
		}

	}
}
void game::locate_pieces() {
	for (int i = 0; i < gridN; ++i) {
		for (int j = 0; j < gridN; ++j) {
			blue_figs[i * 8 + j].setPosition(position_map[i][j].x, position_map[i][j].y);
			pink_figs[i * 8 + j].setPosition(position_map[i][j].x, position_map[i][j].y);
		}
	}
}
void status::change() {
	if (turn == pink) { turn = blue; }
	else if (turn == blue) { turn = pink; }
}
void status::end_check() {
	if (blue_q == 0) { pink_win = 1; }
	else if (pink_q == 0) { blue_win = 1; }
	else if (blue_q + pink_q == 64) {
		if (blue_q > pink_q) { blue_win = 1; }
		else { pink_win = 1; }
	}
}
