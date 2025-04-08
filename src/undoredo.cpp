//-----------------------------------------------------------------------------
// The user-visible undo/redo operation; whenever they change something, we
// record our state and push it on a stack, and we pop the stack when they
// select undo.
//
// Copyright 2008-2013 Jonathan Westhues.
//-----------------------------------------------------------------------------
#include "solvespace.h"

void SolveSpaceUI::UndoRemember() {
    unsaved = true;
    PushFromCurrentOnto(&undo);
    UndoClearStack(&redo);
    UndoEnableMenus();
}

void SolveSpaceUI::UndoUndo() {
    if(undo.cnt <= 0) return;

    PushFromCurrentOnto(&redo);
    PopOntoCurrentFrom(&undo);
    UndoEnableMenus();
}

void SolveSpaceUI::UndoRedo() {
    if(redo.cnt <= 0) return;

    PushFromCurrentOnto(&undo);
    PopOntoCurrentFrom(&redo);
    UndoEnableMenus();
}

void SolveSpaceUI::UndoEnableMenus() {
    SS.GW.undoMenuItem->SetEnabled(undo.cnt > 0);
    SS.GW.redoMenuItem->SetEnabled(redo.cnt > 0);
}

void SolveSpaceUI::PushFromCurrentOnto(UndoStack *uk) {
    if(uk->cnt == MAX_UNDO) {
        UndoClearState(&(uk->d[uk->write]));
        // And then write in to this one again
    } else {
        (uk->cnt)++;
    }

    UndoState *ut = &(uk->d[uk->write]);
    *ut = {};
    ut->group = SK.group;
    for(Group &dest : ut->group) {
        // Zero out all the dynamic stuff that will get regenerated.
        dest.clean = false;
        dest.solved = {};
        dest.polyLoops = {};
        dest.bezierLoops = {};
        dest.bezierOpens = {};
        dest.polyError = {};
        dest.thisMesh = {};
        dest.runningMesh = {};
        dest.thisShell = {};
        dest.runningShell = {};
        dest.displayMesh = {};
        dest.displayOutlines = {};

        dest.impMesh = {};
        dest.impShell = {};
        dest.impEntity = {};
    }
    ut->groupOrder = SK.groupOrder;
    ut->request = SK.request;
    ut->constraint = SK.constraint;
    ut->param = SK.param;
    ut->style = SK.style;
    ut->activeGroup = SS.GW.activeGroup;

    uk->write = WRAP(uk->write + 1, MAX_UNDO);
}

void SolveSpaceUI::PopOntoCurrentFrom(UndoStack *uk) {
    ssassert(uk->cnt > 0, "Cannot pop from empty undo stack");
    (uk->cnt)--;
    uk->write = WRAP(uk->write - 1, MAX_UNDO);

    UndoState *ut = &(uk->d[uk->write]);

    // Do a move of the state from the undo list
    SK.group = std::move(ut->group);
    SK.groupOrder = std::move(ut->groupOrder);
    SK.request = std::move(ut->request);
    SK.constraint = std::move(ut->constraint);
    SK.param = std::move(ut->param);
    SK.style = std::move(ut->style);
    SS.GW.activeGroup = ut->activeGroup;

    *ut = {}; // for good measure

    // And reset the state everywhere else in the program, since the
    // sketch just changed a lot.
    SS.GW.ClearSuper();
    SS.TW.ClearSuper();
    SS.ReloadAllLinked(SS.saveFile);
    SS.GenerateAll(SolveSpaceUI::Generate::ALL);
    SS.ScheduleShowTW();

    // Activate the group that was active before.
    Group *activeGroup = SK.GetGroup(SS.GW.activeGroup);
    activeGroup->Activate();
}

void SolveSpaceUI::UndoClearStack(UndoStack *uk) {
    while(uk->cnt > 0) {
        uk->write = WRAP(uk->write - 1, MAX_UNDO);
        (uk->cnt)--;
        UndoClearState(&(uk->d[uk->write]));
    }
    *uk = {}; // for good measure
}

void SolveSpaceUI::UndoClearState(UndoState *ut) {
    *ut = {};
}

