#pragma once
#include "Loader.h"
#include "board.h"
#include <ctime>
#include<iostream> //for debugging
#include<queue>
#include<vector>

int main() {

	game game0;
	game0.end_window.close();
	while (game0.window.isOpen())
	{
		sf::Event e;
		//°»´ú¨Æ¥ó
		std::queue<sf::Sprite> figs_to_show;
		while (game0.window.pollEvent(e)) {
			if (e.type == sf::Event::Closed) {
				game0.window.close();
				break;
			}
			if (e.type == sf::Event::MouseButtonPressed) {
				//std::cout << e.mouseButton.x << " " << e.mouseButton.y << "\n";
				game0.click_handler(sf::Vector2i(e.mouseButton.x, e.mouseButton.y));
			}
			if (game0.sta.blue_win) {
				game0.end_window.display();
				std::cout << "BLUE WIN~~" << std::endl;
				return 0;
			}
			if (game0.sta.pink_win) {
				game0.end_window.display();
				std::cout << "PINK WIN~~" << std::endl;
				return 0;
			}
		}
		game0.refresh();
		sf::sleep(sf::milliseconds(20));

	}
}