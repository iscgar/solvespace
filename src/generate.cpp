//-----------------------------------------------------------------------------
// Generate our model based on its parametric description, by solving each
// sketch, generating surfaces from the resulting entities, performing any
// requested surface operations (e.g. Booleans) with our model so far, and
// then repeating this process for each subsequent group.
//
// Copyright 2008-2013 Jonathan Westhues.
//-----------------------------------------------------------------------------
#include "solvespace.h"

namespace SolveSpace {

void SolveSpaceUI::MarkGroupDirtyByEntity(hEntity he) {
    Entity *e = SK.GetEntity(he);
    MarkGroupDirty(e->group);
}

void SolveSpaceUI::MarkGroupDirty(hGroup hg, bool onlyThis) {
    bool go = false;
    for(auto const &gh : SK.groupOrder) {
        Group *g = SK.GetGroup(gh);
        if(g->h == hg) {
            go = true;
        }
        if(go) {
            g->clean = false;
            if(onlyThis) break;
        }
    }
    unsaved = true;
    ScheduleGenerateAll();
}

bool Sketch::GroupsInOrder(hGroup hbefore, hGroup hafter) {
    if(hbefore.v != 0 && hafter.v != 0) {
        const Group *before = group.FindByIdNoOops(hbefore);
        if(before == nullptr) {
            return false;
        }
        const Group *after = group.FindByIdNoOops(hafter);
        if(after == nullptr) {
            return false;
        }
        if(before->order >= after->order) {
            return false;
        }
    }
    return true;
}

bool Sketch::EntityExists(hEntity he) {
    // A nonexstient entity is acceptable, though, usually just means it
    // doesn't apply.
    if(he == Entity::NO_ENTITY) return true;
    return entity.FindByIdNoOops(he) ? true : false;
}

bool Sketch::PruneGroup(Group *g) {
    if(GroupsInOrder(g->opA, g->h) &&
       EntityExists(g->predef.origin) &&
       EntityExists(g->predef.entityB) &&
       EntityExists(g->predef.entityC))
    {
        return false;
    }
    group.RemoveById(g->h);
    RegenerateGroupOrder();
    return true;
}


void Sketch::RegenerateGroupOrder() {
    using GroupOrderPair = std::pair<int, hGroup>;
    std::vector<GroupOrderPair> order;
    order.reserve(group.n);
    for(auto &g : group) {
        order.push_back({g.order, g.h});
    }
    std::sort(order.begin(), order.end(), [](const GroupOrderPair &a, GroupOrderPair &b) {
        return a.first < b.first;
    });
    groupOrder.clear();
    groupOrder.reserve(group.n);
    for(auto &p : order) {
        groupOrder.push_back(p.second);
    }
}

Sketch::DeletionStats Sketch::Regenerate(Generate type, bool andFindFree, ParamSet dragged, double chordTol, double *chordTolCalculated) {
    DeletionStats deleted = {};

    RegenerateGroupOrder();
    size_t start = 0, end = groupOrder.size();

    switch(type) {
        case Generate::ALL:
            break;

        case Generate::DIRTY: {
            start = end;

            // Start from the first dirty group, and solve until the active group,
            // since all groups after the active group are hidden.
            // Not using range-for because we're tracking the indices.
            for(size_t i = 0; i < groupOrder.size(); i++) {
                const hGroup hg = groupOrder[i];
                if(start >= end) {
                    const Group *g = GetGroup(hg);
                    if((!g->clean) || !g->IsSolvedOkay()) {
                        start = i;
                    }
                }
                if(hg == SS.GW.activeGroup) {
                    end = i + 1;
                    break;
                }
            }
            if(start < end) {
                SS.nakedEdges.Clear();
            } else {
                // All clean; so just regenerate the entities, and don't solve anything.
            }
            break;
        }

        case Generate::REGEN:
            start = end;
            break;

        case Generate::UNTIL_ACTIVE: {
            end = std::find(groupOrder.begin(), groupOrder.end(), SS.GW.activeGroup) -
                  groupOrder.begin();
            break;
        }
    }

    // Don't lose our numerical guesses when we regenerate.
    ParamList prev = {};
    param.MoveSelfInto(&prev);
    param.ReserveMore(prev.n);
    const int oldEntityCount = entity.n;
    entity.Clear();
    entity.ReserveMore(oldEntityCount);

    // Prune any requests and constraints that belong to non-existent groups
    for(auto &r : request) {
        if(group.FindByIdNoOops(r.group) == nullptr) {
            r.tag = 1;
        }
    }
    for(auto &c : constraint) {
        if(group.FindByIdNoOops(c.group) == nullptr) {
            c.tag = 1;
        }
    }

    // Not using range-for because we're using the index inside the loop.
    for(size_t i = 0; i < groupOrder.size(); ) {
        const hGroup hg = groupOrder[i];
        Group *g = GetGroup(hg);

        // The group may depend on entities or other groups, to define its
        // workplane geometry or for its operands. Those must already exist
        // in a previous group, so check them before generating.
        if(PruneGroup(g)) {
            ++deleted.groups;

            // A group was just removed and the group order was regenerated,
            // so adjust the first and last indices accordingly if we needed
            if(i < end) {
                --end;
            }
            if(i < start) {
                --start;
            }
            // Prune this group's requests
            for(auto &r : request) {
                if(r.group == hg) {
                    r.tag = 1;
                }
            }
            // Prune this group's constraints
            for(auto &c : constraint) {
                if(c.group == hg) {
                    c.tag = 1;
                }
            }
            continue;
        }

        if(hg == Group::HGROUP_REFERENCES) {
            g->varResolutions.clear();
        } else {
            g->varResolutions = GetGroup(groupOrder[i - 1])->varResolutions;
        }

        for (const auto &kv : g->namedParams) {
            g->varResolutions[kv.second.name] = kv.first;
        }

        // Regenerate requests, tagging any dead ones for pruning
        int groupRequestIndex = 0;
        for(auto &r : request) {
            if(r.group != hg) {
                continue;
            }

            // Prune the request if it depends on a pruned workplane
            if(r.workplane != Entity::NO_ENTITY && r.workplane.isFromRequest()) {
                const Request *wr = request.FindByIdNoOops(r.workplane.request());
                if(wr == nullptr || wr->tag) {
                    r.tag = 1;
                    continue;
                }
            }

            r.groupRequestIndex = groupRequestIndex++;
            r.Generate(&entity, &param);
        }

        // Regenerate the group
        g->Generate(&entity, &param);

        // Regenerate constraints, tagging any dead ones for pruning
        for(auto &c : constraint) {
            if(c.group != hg) {
                continue;
            }

            if(!c.IsValid(entity)) {
                if(c.type != Constraint::Type::POINTS_COINCIDENT &&
                   c.type != Constraint::Type::HORIZONTAL &&
                   c.type != Constraint::Type::VERTICAL) {
                    ++deleted.nonTrivialConstraints;
                }

                c.tag = 1;
                continue;
            }

            c.Generate(&param);
        }

        // Use the previous values for params that we've seen before, as
        // initial guesses for the solver.
        for(auto &p : param) {
            Param *newp = &p;
            if(newp->known) {
                continue;
            }

            Param *prevp = prev.FindByIdNoOops(newp->h);
            if(prevp) {
                newp->val = prevp->val;
                newp->free = prevp->free;
            }
        }

        if(hg == Group::HGROUP_REFERENCES) {
            ForceReferences();
            g->solved.how = SolveResult::OKAY;
            g->clean = true;
        } else {
            // this i is an index in groupOrder
            if(i >= start && i < end) {
                // The group falls inside the range, so really solve it,
                // and then regenerate the mesh based on the solved stuff.
                if(chordTolCalculated != nullptr) {
                    SolveGroup(hg, andFindFree, dragged);
                    g->GenerateLoops();
                }
            } else {
                // The group falls outside the range, so just assume that
                // it's good wherever we left it. The mesh is unchanged,
                // and the parameters must be marked as known.
                for(auto &p : param) {
                    Param *newp = &p;

                    Param *prevp = prev.FindByIdNoOops(newp->h);
                    if(prevp) {
                        newp->known = true;
                    }
                }
            }
        }

        ++i;
    }

    prev.Clear();

    // Remove any requests that we tagged for pruning
    const int requests = request.n;
    request.RemoveTagged();
    deleted.requests += requests - request.n;

    // Remove any constraints that we tagged for pruning
    const int constraints = constraint.n;
    constraint.RemoveTagged();
    deleted.constraints += constraints - constraint.n;

    // If we're generating entities for display, first we need to find
    // the bounding box to turn relative chord tolerance to absolute.
    if(chordTolCalculated != nullptr) {
        BBox box = CalculateEntityBBox(/*includeInvisibles=*/true);
        Vector size = box.maxp.Minus(box.minp);
        double maxSize = std::max({ size.x, size.y, size.z });
        *chordTolCalculated = maxSize * chordTol / 100.0;
    }

    // Then generate the shell and mesh
    for(size_t i = start; i < end; ++i) {
        const hGroup hg = groupOrder[i];
        if(hg == Group::HGROUP_REFERENCES) {
            continue;
        }
        Group *g = GetGroup(hg);
        g->GenerateShellAndMesh();
        g->clean = true;
    }

    // And update any reference dimensions with their new values
    for(auto &con : constraint) {
        Constraint *c = &con;
        if(c->reference) {
            c->ModifyToSatisfy(GetGroup(c->group)->varResolutions);
        }
    }

    Platform::FreeAllTemporary();

    return deleted;
}

void Sketch::ForceReferences() {
    // Force the values of the parameters that define the three reference
    // coordinate systems.
    static const struct {
        hRequest    hr;
        Quaternion  q;
    } Quat[] = {
        { Request::HREQUEST_REFERENCE_XY, { 1,    0,    0,    0,   } },
        { Request::HREQUEST_REFERENCE_YZ, { 0.5,  0.5,  0.5,  0.5, } },
        { Request::HREQUEST_REFERENCE_ZX, { 0.5, -0.5, -0.5, -0.5, } },
    };
    for(int i = 0; i < 3; i++) {
        hRequest hr = Quat[i].hr;
        Entity *wrkpl = GetEntity(hr.entity(0));
        // The origin for our coordinate system, always zero
        Entity *origin = GetEntity(wrkpl->point[0]);
        origin->PointForceTo(Vector::From(0, 0, 0));
        origin->construction = true;
        GetParam(origin->param[0])->known = true;
        GetParam(origin->param[1])->known = true;
        GetParam(origin->param[2])->known = true;
        // The quaternion that defines the rotation, from the table.
        Entity *normal = GetEntity(wrkpl->normal);
        normal->NormalForceTo(Quat[i].q);
        GetParam(normal->param[0])->known = true;
        GetParam(normal->param[1])->known = true;
        GetParam(normal->param[2])->known = true;
        GetParam(normal->param[3])->known = true;
    }
}

void SolveSpaceUI::GenerateAll(Generate type, bool andFindFree) {
    const uint64_t startMillis = GetMilliseconds();

    int first = 0, last = 0, i;

    const auto deleted = SK.Regenerate(type, andFindFree, GetDraggedParams(), chordTol,
                                       !exportMode ? &chordTolCalculated : nullptr);

    for(hGroup hg : SK.groupOrder) {
        Group *g = SK.GetGroup(hg);
        if(!g->IsSolvedOkay()) {
            TextWindow::ReportHowGroupSolved(g->h);
        }
    }

    // Make sure the point that we're tracing exists.
    if(traced.point.v) {
        if(!SK.entity.FindByIdNoOops(traced.point)) {
            traced.point = Entity::NO_ENTITY;
        } else {
            // And if we're tracing a point, add its new value to the path
            Entity *pt = SK.GetEntity(traced.point);
            traced.path.AddPoint(pt->PointGetNum());
        }
    }

    GW.Invalidate();

    // Remove nonexistent selection items, for same reason we waited till
    // the end to put up a dialog box.
    GW.ClearNonexistentSelectionItems();

    const uint64_t endMillis = GetMilliseconds();

    if(endMillis - startMillis > 30) {
        const char *typeStr = "";
        switch(type) {
            case Generate::DIRTY:           typeStr = "DIRTY";        break;
            case Generate::ALL:             typeStr = "ALL";          break;
            case Generate::REGEN:           typeStr = "REGEN";        break;
            case Generate::UNTIL_ACTIVE:    typeStr = "UNTIL_ACTIVE"; break;
        }
        if(endMillis)
        dbp("Generate::%s took %lld ms",
            typeStr,
            GetMilliseconds() - startMillis);
    }

    if(deleted.requests > 0 || deleted.constraints > 0 || deleted.groups > 0) {
        // All sorts of interesting things could have happened; for example,
        // the active group or active workplane could have been deleted. So
        // clear all that out.
        if(deleted.groups > 0) {
            SS.TW.ClearSuper();
        }
        ScheduleShowTW();
        GW.ClearSuper();

        // People get annoyed if I complain whenever they delete any request,
        // and I otherwise will, since those always come with pt-coincident
        // constraints.
        if(deleted.requests > 0 || deleted.nonTrivialConstraints > 0 ||
           deleted.groups > 0)
        {
            // Don't display any errors until we've regenerated fully. The
            // sketch is not necessarily in a consistent state until we've
            // pruned any orphaned etc. objects, and the message loop for the
            // messagebox could allow us to repaint and crash. But now we must
            // be fine.
            Message("Additional sketch elements were deleted, because they "
                    "depend on the element that was just deleted explicitly. "
                    "These include: \n"
                    "     %d request%s\n"
                    "     %d constraint%s\n"
                    "     %d group%s"
                    "%s",
                       deleted.requests, deleted.requests == 1 ? "" : "s",
                       deleted.constraints, deleted.constraints == 1 ? "" : "s",
                       deleted.groups, deleted.groups == 1 ? "" : "s",
                       undo.cnt > 0 ? "\n\nChoose Edit -> Undo to undelete all elements." : "");
        }
    }

    allConsistent = true;
    SS.GW.persistentDirty = true;
    SS.centerOfMass.dirty = true;
}

void SolveSpaceUI::UpdateCenterOfMass() {
    SMesh *m = &(SK.GetGroup(SS.GW.activeGroup)->displayMesh);
    SS.centerOfMass.position = m->GetCenterOfMass();
    SS.centerOfMass.dirty = false;
}

ParamSet SolveSpaceUI::GetDraggedParams() const {
    ParamSet dragged;

    for(int i = -1; i < SS.GW.pending.points.n; i++) {
        hEntity hp;
        if(i == -1) {
            hp = SS.GW.pending.point;
        } else {
            hp = SS.GW.pending.points[i];
        }
        if(!hp.v) continue;

        // The pending point could be one in a group that has not yet
        // been processed, in which case the lookup will fail; but
        // that's not an error.
        Entity *pt = SK.entity.FindByIdNoOops(hp);
        if(pt) {
            switch(pt->type) {
                case Entity::Type::POINT_N_TRANS:
                case Entity::Type::POINT_IN_3D:
                case Entity::Type::POINT_N_ROT_AXIS_TRANS:
                    dragged.insert(pt->param[0]);
                    dragged.insert(pt->param[1]);
                    dragged.insert(pt->param[2]);
                    break;

                case Entity::Type::POINT_IN_2D:
                    dragged.insert(pt->param[0]);
                    dragged.insert(pt->param[1]);
                    break;

                default: // Only the entities above can be dragged.
                    break;
            }
        }
    }
    if(SS.GW.pending.circle.v) {
        Entity *circ = SK.entity.FindByIdNoOops(SS.GW.pending.circle);
        if(circ) {
            Entity *dist = SK.GetEntity(circ->distance);
            switch(dist->type) {
                case Entity::Type::DISTANCE:
                    dragged.insert(dist->param[0]);
                    break;

                default: // Only the entities above can be dragged.
                    break;
            }
        }
    }
    if(SS.GW.pending.normal.v) {
        Entity *norm = SK.entity.FindByIdNoOops(SS.GW.pending.normal);
        if(norm) {
            switch(norm->type) {
                case Entity::Type::NORMAL_IN_3D:
                    dragged.insert(norm->param[0]);
                    dragged.insert(norm->param[1]);
                    dragged.insert(norm->param[2]);
                    dragged.insert(norm->param[3]);
                    break;

                default: // Only the entities above can be dragged.
                    break;
            }
        }
    }

    return dragged;
}

void Sketch::WriteEqSystemForGroup(hGroup hg, ParamSet dragged) {
    // Clear out the system to be solved.
    sys.entity.Clear();
    sys.param.Clear();
    sys.eq.Clear();
    // And generate all the params for requests in this group
    for(auto &req : request) {
        Request *r = &req;
        if(r->group != hg || r->tag) continue;

        r->Generate(&(sys.entity), &(sys.param));
    }
    for(auto &con : constraint) {
        Constraint *c = &con;
        if(c->group != hg || c->tag) continue;

        c->Generate(&(sys.param));
    }
    // And for the group itself
    Group *g = GetGroup(hg);
    g->Generate(&(sys.entity), &(sys.param));
    // Set the initial guesses for all the params
    for(auto &param : sys.param) {
        Param *p = &param;
        const auto namedParam = g->GetNamedParam(p->h);
        if(namedParam.first == p->h && namedParam.second.locked) {
            p->known = true;
        } else {
            p->known = false;
        }
        p->val = GetParam(p->h)->val;
    }

    sys.dragged = std::move(dragged);
}

void Sketch::SolveGroup(hGroup hg, bool andFindFree, ParamSet dragged) {
    WriteEqSystemForGroup(hg, std::move(dragged));
    Group *g = GetGroup(hg);
    g->solved.remove.clear();
    g->solved.findToFixTimeout = SS.timeoutRedundantConstr;
    SolveResult how = sys.Solve(g, &(g->solved.dof),
                                   &(g->solved.remove),
                                   /*andFindBad=*/!g->allowRedundant,
                                   /*andFindFree=*/andFindFree,
                                   /*forceDofCheck=*/!g->dofCheckOk);
    if(how == SolveResult::OKAY) {
        g->dofCheckOk = true;
    }
    g->solved.how = how;
    Platform::FreeAllTemporary();
}

SolveResult Sketch::TestRankForGroup(hGroup hg, ParamSet dragged, int *rank) {
    Group *g = GetGroup(hg);
    // If we don't calculate dof or redundant is allowed, there is
    // no point to solve rank because this result is not meaningful
    if(g->suppressDofCalculation || g->allowRedundant) return SolveResult::OKAY;
    WriteEqSystemForGroup(hg, std::move(dragged));
    SolveResult result = sys.SolveRank(g, rank);
    Platform::FreeAllTemporary();
    return result;
}

bool SolveSpaceUI::ActiveGroupsOkay() {
    for(hGroup hg : SK.groupOrder) {
        if(!SK.GetGroup(hg)->IsSolvedOkay()) {
            return false;
        }
        if(hg == GW.activeGroup) {
            break;
        }
    }
    return true;
}

} // namespace SolveSpace
