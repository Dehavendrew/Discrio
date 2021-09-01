// Discrio.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#define OLC_PGE_APPLICATION
#include <iostream>
#include "C:\Users\Drew\Desktop\olcPixelGameEngine-master\olcPixelGameEngine.h"
#include "C:\Users\Drew\Desktop\synth-master\olcNoiseMaker.h"
#include <vector>
#include <cmath>
#include <stdio.h>   
#include <stdlib.h>     
#include <time.h>  
#include <cstdlib>
#include <atomic>
#include <condition_variable>
#include <random>

constexpr int nMaxThreads = 8;

int nMode = 0;

float fPrevTime = 0.0f;
float fMasterTime = 0.0f;

double dAngFreq = 8;
double dCAngFreq = 2;
double dWaveNum = .03;

double dAmplitude = 0;
double dDCOffset = 0;
double dMiddleScreen = 0;
int nAmplitudeScale = 8;
int nOffsetScale = 2;

double dTimeDecay = 0;//-0.01;
double dSpaceDecay = 0;//-0.001;
double dSpaceFreqMod = 0; // .009;
double dTimeFreqMod = 0; //0.9;
double dRadiusFreq = 0;
double dCircphase = 0;
bool bSyncTimeFreq = true;

int nHeight = 0;
int nWidth = 0;
int nHarmonicOrder = 1;
int nShape = 2;
int nRadius = 20;

double dRedPhase = 0;
double dGreenPhase = 0;
double dBluePhase = 0;

double dNoise = 0.0;
vector<double> dPrevNoise;


std::atomic<double> dSyncTimeFreq = 0;
//double dSyncTimeFerq = 0;
std::mutex m;

std::random_device mr;
std::mt19937 mmt(mr());
std::uniform_real_distribution<double> dist(1.0, 1.0);

void randomize() {

	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(0.0, 1.0);
	std::uniform_real_distribution<double> dist1(2.0, 10.0);
	std::uniform_real_distribution<double> dist2(0.0, 0.003);
	std::uniform_real_distribution<double> dist3(0.0, 6.0);
	std::uniform_real_distribution<double> dist4(0.0, 0.5);
	std::uniform_real_distribution<double> dist5(0.0, 0.001);

	dAngFreq = dist1(mt);
	dCAngFreq = dist1(mt);
	dWaveNum = dist4(mt);

	nAmplitudeScale = dist1(mt);
	dAmplitude = nHeight / 3 * dist(mt);

	dTimeDecay = -1 * dist2(mt);//-0.01;
	dSpaceDecay = dist2(mt);//-0.001;
	if (dist3(mt) > 2.5) {
		dSpaceDecay *= -1;
	}
	else dSpaceDecay /= 2;

	dSpaceFreqMod = dist5(mt); // .009;
	dTimeFreqMod = dist(mt); //0.9;
	dRadiusFreq = dist(mt);

	if (dist3(mt) > 3.0) {
		bSyncTimeFreq = true;
	}
	else bSyncTimeFreq = false;

	nHarmonicOrder = dist1(mt);
	nShape = dist1(mt);

	dRedPhase = dist3(mt);
	dGreenPhase = dist3(mt);
	dBluePhase = dist3(mt);

}


class Example : public olc::PixelGameEngine
{
public:

	Example()
	{
		sAppName = "Example";
	}


	std::vector<std::vector<int>> points;


public:

	bool OnUserCreate() override
	{
		// Called once at the start, so create things here
		nHeight = ScreenHeight();
		nWidth = ScreenWidth();
		dAmplitude = ScreenHeight() / nAmplitudeScale;
		dDCOffset = ScreenHeight() / nOffsetScale;
		dMiddleScreen = ScreenWidth() / 2;

		if(nMode == 2) tree(ScreenWidth(),ScreenHeight());

		for (int i = 0; i < nHarmonicOrder; ++i) {
			dPrevNoise.push_back(1.0);
		}
		
		return true;
	}

	void tree(int w, int h) {

		std::vector<int> start;
		start.push_back(ScreenHeight() / 2);
		points.push_back(start);

		int nProbFork = 210;

		for (int i = 0; i <= w; ++i) {
			std::vector<int> row;
			for (int j = 0; j < points[i].size(); ++j) {

				int chance = rand() % nProbFork + 1;
				if (chance == 100) {
					row.push_back(points[i][j] + 1);
					row.push_back(points[i][j] - 1);
				}
				else if (chance % 3 == 0) row.push_back(points[i][j] + 1);
				else if (chance % 3 == 1) row.push_back(points[i][j] - 1);
				else if (chance % 3 == 2) row.push_back(points[i][j]);

				if (row.size() >= h) continue;
			}

			points.push_back(row);

		}

		for (int i = 0; i < points.size(); ++i) {
			for (int j = 0; j < points[i].size(); ++j) {
				Draw(i, points[i][j], olc::Pixel(255,255,0));
			}
		}

	}

	double wave(float x, float t, vector<double> noise) {

		double dVal = dDCOffset;

	//	m.lock();
		dSyncTimeFreq.store(dTimeFreqMod * t, std::memory_order_release);
		//dSyncTimeFreq = dTimeFreqMod * t;
		//m.unlock();

		for (int i = 1; i <= nHarmonicOrder; ++i) {
			
			if (i % nShape == 0) continue;
			
			int ntimeFreq = 1;

			if (bSyncTimeFreq) {
				ntimeFreq = i;
			}

		//	dVal += exp(dTimeDecay * t + dSpaceDecay * x) * (dAmplitude * (1.0 / i)* sin(cos(dSyncTimeFreq + dSpaceFreqMod * x) * dWaveNum * i * x + dAngFreq * ntimeFreq * t * noise[i - 1]));
			dVal += exp(dTimeDecay * t + dSpaceDecay * x) * (dAmplitude * (1.0 / i)* sin(cos(dSyncTimeFreq + dSpaceFreqMod * x) * dWaveNum * i * x) * sin(dAngFreq * ntimeFreq * t * noise[i - 1]));


		}


		return dVal;

	}

	olc::Pixel waveColor(float t) {

		int color [3];

		color[0] = abs(255 * sin(dCAngFreq * t + dRedPhase));
		color[1] = abs(255 * sin(dCAngFreq * t + dGreenPhase));
		color[2] = abs(255 * sin(dCAngFreq * t + dBluePhase));

		olc::Pixel pixel(color[0], color[1], color[2]);

		return pixel;
	}

	void drawSeg(int x, float t, vector<double> noise, bool clear=false) {

		int d1;
		int d2;

		
		d1 = wave(x, t, noise);
		d2 = wave(x + 1, t, noise);


		olc::Pixel pixel;

		int* colors;
		if (clear) {
			pixel = olc::Pixel(0,0,0);
		}
		else {
			pixel = waveColor(t);
		}

		if (d1 < d2) {
			for (int i = d1; i <= d2; ++i) {

				Draw(x, i, pixel);

			}
		}

		else {
			for (int i = d2; i <= d1; ++i) {

				Draw(x, i, pixel);

			}
		}
		
	}

	void drawCirc(float x, float t, vector<double> noise, bool clear = false) {


		int d1;
		int d2;
		

		d1 = wave(x, t, noise);
		d2 = wave(x + 1, t, noise);

		olc::Pixel pixel;

		int* colors;
		if (clear) {
			pixel = olc::Pixel(0, 0, 0);
		}
		else { 
			pixel = waveColor(t);
		}

		float radius1 = nRadius * cos(t * dRadiusFreq) + d1 - dDCOffset + 10;
		float radius2 = nRadius * cos(t * dRadiusFreq) + d2 - dDCOffset + 10;

		int y1 = radius1 * sin(x / ScreenWidth() * 2 * 3.14159 + dCircphase) + dDCOffset;
		int y2 = radius2 * sin((x + 1) / ScreenWidth() * 2 * 3.14159 + dCircphase) + dDCOffset;

		int x1 = radius1 * cos(x / ScreenWidth() * 2 * 3.14159 + dCircphase) + dMiddleScreen;
		int x2 = radius2 * cos((x + 1) / ScreenWidth() * 2 * 3.14159 + dCircphase) + dMiddleScreen;
		
		DrawLine(olc::vi2d(x1, y1), olc::vi2d(x2, y2), pixel);

	}

	olc::Pixel eval(int x, int y, int t, int mode) {

		if (mode == 0) {
			//Grid Mode
			if (x % 50 == 0 || y % 50 == 0) {
				return olc::Pixel(abs(x) % 256, abs(y) % 256, 255);
			}
		}
		if (mode == 1) {
			//Sine
			double dSpaceFreq = 0.001 * t;

			int r = abs(255 * sin(dSpaceFreq * (pow(x, 2) + pow(y, 2)) + 0));
			int g = abs(255 * sin(dSpaceFreq * (pow(x, 2) + pow(y, 2)) + 2));
			int b = abs(255 * sin(dSpaceFreq * (pow(x, 2) + pow(y, 2)) + 4));
			return olc::Pixel(r, g, b);
		}

		return olc::Pixel(0, 0, 0);
	}

	int* transform(int x, int y, int type = 0) {

		int points [2];

		points[0] = x - ScreenWidth() / 2;
		points[1] = y - ScreenHeight() / 2;

		if (type == 1) {
			int a = 1;
			int b = 0;
			int c = 1;
			int d = 1;

			points[0] = a * points[0] + b * points[1];
			points[1] = c * points[0] + d * points[1];
		}

		return points;

	}

	void paintSection(int start, int stop) {
		for (int i = start; i < stop; ++i) {
			for (int j = 0; j < ScreenHeight(); ++j) {

				int* points = transform(i, j);

				olc::Pixel color = eval(points[0], points[1], fMasterTime, 1);

				Draw(i, j, color);
			}
		}
	}

	void unrandomize() {
		Clear(olc::BLACK);
		dAngFreq = 1.1;
		dCAngFreq = 0.001;
		dWaveNum = 0.03;

		dAmplitude = ScreenHeight() / 7;

		dTimeDecay = -0.00001;
		dSpaceDecay = -0.00001;
		dSpaceFreqMod = 0.00001;
		dTimeFreqMod = 0.00001;
		dRadiusFreq = 0.00001;
		dCircphase = 0;
		nRadius = 50;

		
		bSyncTimeFreq = true;

		nHarmonicOrder = 1;
		nShape = 5;

		dRedPhase = 2;
		dGreenPhase = 2;
		dBluePhase = 2;

		dNoise = 0;
		std::uniform_real_distribution<double> dist1(1.0, 1.0);
		dist = dist1;
	}


	bool OnUserUpdate(float fElapsedTime) override
	{
		if (GetKey(olc::Key::ENTER).bPressed) {
			randomize();
			while (dPrevNoise.size() < nHarmonicOrder) {
				dPrevNoise.push_back(1.0);
			}
		}
		if (GetKey(olc::Key::SHIFT).bHeld) {
			Clear(olc::BLACK);
		}
		if (GetKey(olc::Key::Q).bPressed) {
			bSyncTimeFreq = bSyncTimeFreq ^ true;
		}
		if (GetKey(olc::Key::K1).bHeld) {
			nMode = 0;
		}
		if (GetKey(olc::Key::K2).bHeld) {
			nMode = 1;
		}
		if (GetKey(olc::Key::K3).bHeld) {
			nMode = 3;
		}
		if (GetKey(olc::Key::K4).bHeld) {
			nMode = 4;
		}
		if (GetKey(olc::Key::K5).bHeld) {
			nMode = 5;
		}
		if (GetKey(olc::Key::K6).bHeld) {
			nMode = 6;
		}

		if (GetKey(olc::Key::BACK).bHeld) {
			unrandomize();
		}

		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::A).bHeld) {
			dAmplitude *= 1.01;
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::A).bHeld) {
			dAmplitude /= 1.01;
		}

		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::X).bHeld) {
			dSpaceDecay -= .00001;
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::X).bHeld) {
			dSpaceDecay += .00001;
		}

		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::M).bHeld) {
			dSpaceFreqMod *= 1.01;
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::M).bHeld) {
			dSpaceFreqMod /= 1.01;
		}

		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::U).bHeld) {
			dTimeFreqMod *= 1.01;
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::U).bHeld) {
			dTimeFreqMod /= 1.01;
		}
		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::T).bHeld) {
			dRadiusFreq *= 1.01;
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::T).bHeld) {
			dRadiusFreq /= 1.01;
		}
		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::E).bHeld) {
			dCircphase += .01;
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::E).bHeld) {
			dCircphase -= .01;
		}

		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::K).bHeld) {
			dWaveNum *= 1.01;
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::K).bHeld) {
			dWaveNum /= 1.01;
		}

		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::F).bHeld) {
			dAngFreq *= 1.01;
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::F).bHeld) {
			dAngFreq /= 1.01;
		}

		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::C).bHeld) {
			dCAngFreq *= 1.01;
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::C).bHeld) {
			dCAngFreq /= 1.01;
		}

		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::N).bPressed) {
			++nHarmonicOrder;
			dPrevNoise.push_back(1.0);
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::N).bPressed) {
			--nHarmonicOrder;
		}

		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::R).bHeld) {
			++nRadius;
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::R).bHeld) {
			if(nRadius != 1) --nRadius;
		}

		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::S).bPressed) {
			++nShape;
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::S).bPressed) {
			if(nShape != 2) --nShape;
		}

		if (GetKey(olc::Key::I).bHeld) {
			dRedPhase += 0.1;
		}
		if (GetKey(olc::Key::O).bHeld) {
			dGreenPhase += 0.1;
		}
		if (GetKey(olc::Key::P).bHeld) {
			dBluePhase += 0.1;
		}

		if (GetKey(olc::Key::UP).bHeld && GetKey(olc::Key::W).bHeld) {
			dNoise += 0.0001;
			std::uniform_real_distribution<double> dist1(1.0 - dNoise, 1.0 + dNoise);
			dist = dist1;
		}
		if (GetKey(olc::Key::DOWN).bHeld && GetKey(olc::Key::W).bHeld) {
			std::uniform_real_distribution<double> dist1(1.0 - dNoise, 1.0 + dNoise);
			dist = dist1;
			dNoise -= 0.0001;
			if (dNoise < 0) dNoise = 0;
		}
		if (GetKey(olc::Key::TAB).bPressed) {
			std::cout << "Synced Wave: " << bSyncTimeFreq << endl;
			std::cout << "nMode: " << nMode << endl;
			std::cout << "dAmplitude: " << dAmplitude << endl;
			std::cout << "dSpaceDecay: " << dSpaceDecay << endl;
			std::cout << "dSpaceFreqMod: " << dSpaceFreqMod << endl;
			std::cout << "dTimeFreqMod: " << dTimeFreqMod << endl;
			std::cout << "dWaveNum: " << dWaveNum << endl;
			std::cout << "dAngFreq: " << dAngFreq << endl;
			std::cout << "dCAngFreq: " << dCAngFreq << endl;
			std::cout << "nHarmonicOrder: " << nHarmonicOrder << endl;
			std::cout << "nShape: " << nShape << endl;
			std::cout << "dRedPhase: " << dRedPhase << endl;
			std::cout << "dGreenPhase: " << dGreenPhase << endl;
			std::cout << "dBluePhase: " << dBluePhase << endl;
			std::cout << "dNoise: " << dNoise << endl;
			std::cout << endl << endl << endl;
		}


		vector<double> dCurrNoise;

		for (int k = 0; k < nHarmonicOrder; ++k) {
			dCurrNoise.push_back(dist(mmt));
		}

		fMasterTime += fElapsedTime;

		if (nMode == 0 || nMode == 1) {
			for (int x = 0; x < ScreenWidth(); ++x) {

				if(nMode == 0) drawSeg(x, fPrevTime, dPrevNoise, 1);

				drawSeg(x, fMasterTime, dCurrNoise);

			}
		}

		if (nMode == 5 || nMode == 6) {
			for (int x = 0; x < ScreenWidth(); ++x) {

				if (nMode == 5)drawCirc(x, fPrevTime, dCurrNoise, 1); 

				drawCirc(x, fMasterTime, dCurrNoise);

			}
		}

		if (nMode == 3) {

			int nSectionWidth = ScreenWidth() / nMaxThreads;
			std::thread t[nMaxThreads];

			for (int i = 0; i < nMaxThreads; ++i) {
				t[i] = std::thread(&Example::paintSection, this, nSectionWidth * i, nSectionWidth * (i + 1));
			}

			for (int i = 0; i < nMaxThreads; i++)
				t[i].join();

		}

		if (nMode == 4) {

			for (int i = 0; i < ScreenWidth(); ++i) {
				for (int j = 0; j < ScreenHeight(); ++j) {

					int* points = transform(i, j);

					olc::Pixel color = eval(points[0], points[1], fMasterTime, 1);

					Draw(i, j, color);
				}
			}

		}

		fPrevTime = fMasterTime;
		dPrevNoise = dCurrNoise;

		return true;
	}


};

double MakeNoise(double dTime)
{
	double dVal = 0;

	for (int i = 1; i <= nHarmonicOrder; ++i) {

		if (i % nShape == 0) continue;

		int ntimeFreq = 1;

		if (bSyncTimeFreq) {
			ntimeFreq = i;
		}

		double noise = 1.0;
		if (dPrevNoise.size() > i) {
			noise = dPrevNoise[i];
		}

		dVal += exp(dTimeDecay * dTime) * ((1.0 / i) * sin(cos(dTimeFreqMod * dTime/*dSyncTimeFreq.load(std::memory_order_relaxed)*/) * dAngFreq * 25 * 2.0 * 3.14159 * ntimeFreq * dTime * noise));

	}

	return dVal * 0.2 * (dAmplitude / nHeight); // Master Volume
}

void sound() {
	vector<wstring> devices = olcNoiseMaker<short>::Enumerate();
	olcNoiseMaker<short> sound(devices[0], 44100, 1, 8, 512);
	sound.SetUserFunction(MakeNoise);

	while (1) {

	}
}

void show() {
	randomize();
	Example demo;

	if (demo.Construct(1200, 600, 1, 1))
		demo.Start();
}

int main()
{
	bool sounds = false;
	thread threads[2];

	if (sounds) {
		threads[0] = std::thread(sound);
		threads[1] = std::thread(show);

		threads[1].join();
		threads[0].join();
	}

	
	else show();

	return 0;
}

