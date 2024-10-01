#include "pch.h"
#include "utils.h"
#include "../msmodel/event_loop.h"
using namespace std;


/* TODO :
pass this function into EventLoop member function (in event_loop.cpp - projet msmodel)
once completely tested */
void EventLoop::randomStep() {

	typedef void(EventLoop::* fptr)();
	std::vector<fptr> fvec = { &EventLoop::cellBirthAll, &EventLoop::cellDeathAll,
							   &EventLoop::plasmidTransferAll, &EventLoop::plasmidLossAll,
							   &EventLoop::cellDiffusionAll };
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(fvec.begin(), fvec.end(), g);
	for (int i = 0; i < fvec.size(); i++) {
		(this->*fvec[i])();
	};

	u.shuffle();

}


TEST(EventLoop, randomStep) {

	Universe u;
	setup_01(u);
	Sampler sp;
	EventLoop ev(u, sp);

	ev.randomStep();

}