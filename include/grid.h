/*
    Grid class
*/
#ifndef GRID_H
#define GRID_H

class Grid
{
private:
    /* Flag to indicate a 2d grid or 3d grid */
    const bool is_2d_;

    /* Number of grids on x, y, and z axis */
    const unsigned X_;
    const unsigned Y_;
    const unsigned Z_;

    /* Size of each grid */
    const unsigned INTERVAL_X_;
    const unsigned INTERVAL_Y_;
    const unsigned INTERVAL_Z_;

public:
    /* Default constructor */
    Grid() : is_2d_(true), X_(0), Y_(0), Z_(0), INTERVAL_X_(0), INTERVAL_Y_(0), INTERVAL_Z_(0) {}

    Grid(const bool is_2d, const unsigned X, const unsigned Y, const unsigned Z, const unsigned INTERVAL_X, const unsigned INTERVAL_Y, const unsigned INTERVAL_Z);

    /* Destructor to clean up allocated data structures */
    ~Grid()
    {
    }
};

#endif