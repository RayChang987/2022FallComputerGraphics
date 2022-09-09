#pragma once;
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include "board.h"


//every object on the chess board
class status {
public:
	int blue_q = 2, pink_q = 2;
	const int blue = 1, pink = 2, empty = 0;
	int turn = blue;
	bool blue_win = false, pink_win = false;
	void change();
	void end_check();
};
class game {
public:
	game() {
		main_board.board_init();
		
		//create the window 
		window.create(sf::VideoMode(winW, winH), "Reversi2209");
		end_window.create(sf::VideoMode(100, 100), "END!!!");

		//load materials
		
		buffer.loadFromFile("Sound/chess-move.wav");
		background_text.loadFromFile("Images/board.png");
		blue_text.loadFromFile("Images/blue.png");
		pink_text.loadFromFile("Images/pink.png");

		//build figures
		background = sf::Sprite(background_text);
		for(int i = 0; i<64; ++i){
			blue_figs[i].setTexture(blue_text);
			pink_figs[i].setTexture(pink_text);
		}
		//sound
		sound.setBuffer(buffer);

		//calculate the position of pieces
		fill_position_map();

		//puts all pieces on the board
		locate_pieces();

	}
	// App object
	//Sound
	sf::SoundBuffer buffer;
	sf::Sound sound;

	//Graphics
	void read_board_and_react();
	sf::RenderWindow window;
	sf::RenderWindow end_window;
	static const int winW = 560, winH = 560, gridN = 8;
	sf::Texture background_text, blue_text, pink_text;
	sf::Sprite background;
	sf::Sprite blue_figs[64], pink_figs[64];
	
	//game
	board main_board;
	void click_handler(const sf::Vector2i&);
	void fill_position_map();
	void locate_pieces();
	void refresh();
	status sta;

	sf::Vector2i position_map[gridN][gridN];
};



