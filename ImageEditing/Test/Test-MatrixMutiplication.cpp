#include<iostream>
using namespace std;

int main(){
    int r = (103>>3), g = (113>>3), b = (19>>3);
    int encoded_data = ((r <<10)) + ((g <<5)) + ((b<<0));
    int rr = encoded_data >> 10;
    encoded_data -= (encoded_data>>10) << 10;
    int gg = encoded_data >> 5;
    encoded_data -= (gg<<5);
    int bb= encoded_data;
    cout<< r<<" "<<g<<" "<<b<<endl;
    cout<< rr<<" "<<gg<<" "<<bb<<endl;
}
