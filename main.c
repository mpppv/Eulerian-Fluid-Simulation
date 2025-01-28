//
//  main.c
//  Fluid_Simulation
//
//  Created by mpppv on 2025/01/26.
//
#include <SDL2/SDL.h>
#include <stdio.h>
#include <stdbool.h>

#define COLOR_BLACK         0x00000000
#define COLOR_GREY          0x0f0f0f0f
#define LIQUID_COLOR        0x34c3eb
#define SOLID_COLOR         0xbea1bf

#define PRESSURE_MODE       true

#define WINDOW_WIDTH        900
#define WINDOW_HEIGHT       600

#define SPHERE_RADIUS       30

#define CELL_SIZE           5
#define LINE_WIDTH          2
#define COLUMN_NUM          WINDOW_WIDTH / CELL_SIZE
#define ROW_NUM             WINDOW_HEIGHT / CELL_SIZE
#define GRID_COLOR          0x262626

#define CELL_SIZE_M         0.3
#define TIME_STEP_HZ        30
#define GRAVITY             9.81
#define DENSITY             1000
#define OVERRELAXATION      1.7


enum matter_type{
    liquid,
    solid
};

typedef struct{
    enum matter_type type;
    int x, y;
    
    double v, u;
    double ro;
    double p;
}cell;

typedef struct{
    int x, y;
    double v, u;
    
    bool    ocup_cll[ROW_NUM][COLUMN_NUM];
}solid_obj;

cell** init_cells(void);
solid_obj init_obj(void);
void set_cell_to(cell* cell, enum matter_type type);
void update_obj(solid_obj* obj, int x, int y);

void add_gravity(SDL_Surface* surface, cell** cells);
void projection(SDL_Surface* surface, cell** cells);
void advection(SDL_Surface* surface, cell** cells);

void draw_grid(SDL_Surface* surface);
void draw_cell(SDL_Surface* surface, int x, int y, Uint32 color);
void draw_cells(SDL_Surface* surface, cell** cells);

int main(void) {
    
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Init(SDL_INIT_EVENTS);
    SDL_Window* window = SDL_CreateWindow("Liquid Simulation",
                                          SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                          WINDOW_WIDTH,           WINDOW_HEIGHT, 0);
    SDL_Surface* surface = SDL_GetWindowSurface(window);
    
    solid_obj obj = init_obj();
    cell** cells = init_cells();
    
    
    bool running = true;
    int mouse_x = 0, mouse_y = 0;
    SDL_Event event;
    while(running){
        while(SDL_PollEvent(&event)){
            if(event.type == SDL_QUIT) running = false;
            else if(event.button.button == SDL_BUTTON_LEFT){
                SDL_GetMouseState(&mouse_x, &mouse_y);
                update_obj(&obj, mouse_x/CELL_SIZE, ROW_NUM-1 - mouse_y/CELL_SIZE);
                //SDL_GetModState() & KMOD_SHIFT
            }
        }
        
        
        SDL_Rect black_rect = (SDL_Rect){0, 0, WINDOW_WIDTH, WINDOW_HEIGHT};
        SDL_FillRect(surface, &black_rect, COLOR_BLACK);
        draw_cells(surface, cells);
        //draw_grid(surface);
        
        add_gravity(surface, cells);
        projection(surface, cells);
        advection(surface, cells);
        
//        for(int i = 0; i < ROW_NUM; ++i){
//            for(int j = 0; j < COLUMN_NUM; ++j){
//                printf("%---d\t", (int)(cells[ROW_NUM-1-i][j].p/10000));
//            }
//            printf("\n");
//        }
//        printf("--------------------------------------------------------\n");
        
        
        SDL_UpdateWindowSurface(window);
        SDL_Delay(1000/TIME_STEP_HZ);
    }
}


cell** init_cells(void){
    cell** cells = (cell**)malloc(ROW_NUM * sizeof(cell*));
    
    for (int i = 0; i < ROW_NUM; ++i) {
        cells[i] = (cell*)malloc(COLUMN_NUM * sizeof(cell));
        for (int j = 0; j < COLUMN_NUM; ++j)
            cells[i][j] = (cell){liquid, i, j, 0, 0, 0, 0};
    }
    return cells;
}

solid_obj init_obj(void){
    bool** posit = (bool**)malloc(ROW_NUM * sizeof(bool*));
    
    for (int i = 0; i < ROW_NUM; ++i) {
        posit[i] = (bool*)malloc(COLUMN_NUM * sizeof(bool));
        for (int j = 0; j < COLUMN_NUM; ++j)
            posit[i][j] = (i - ROW_NUM/2)*(i - ROW_NUM/2) + (j - COLUMN_NUM/2)*(j - COLUMN_NUM/2) <= SPHERE_RADIUS * SPHERE_RADIUS;
    }
    
    return (solid_obj){ROW_NUM/2, COLUMN_NUM/2, 0, 0, posit};
}

void set_cell_to(cell* cell, enum matter_type type){
    cell->type = type;
    cell->v = cell->u = 0;
    cell->ro = cell->p = 0;
}

void update_obj(solid_obj* obj, int x, int y){
    
}

void add_gravity(SDL_Surface* surface, cell** cells){
    for(int i = 0; i < ROW_NUM; ++i)
        for(int j = 0; j < COLUMN_NUM; ++j)
            if(cells[i][j].type == liquid)
                cells[i][j].v -= GRAVITY / TIME_STEP_HZ;
}

void projection(SDL_Surface* surface, cell** cells){
    double d;
    int s, s1, s2, s3, s4;
    for(int i = 0; i < ROW_NUM; ++i)
        for(int j = 0; j < COLUMN_NUM; ++j)
            if(cells[i][j].type == liquid){
                if(i == ROW_NUM-1 && j == COLUMN_NUM-1)
                    d = cells[i][j].u - cells[i][j].v;
                else if(i == ROW_NUM-1)
                    d = cells[i][j+1].u - cells[i][j].u - cells[i][j].v;
                else if(j == COLUMN_NUM-1)
                    d = -cells[i][j].u + cells[i+1][j].v - cells[i][j].v;
                else
                    d = cells[i][j+1].u - cells[i][j].u + cells[i+1][j].v - cells[i][j].v;
                d *= OVERRELAXATION;
                
                s1 = (i > 0)            && (cells[i-1][j].type == liquid);
                s2 = (i < ROW_NUM-1)    && (cells[i+1][j].type == liquid);
                s3 = (j > 0)            && (cells[i][j-1].type == liquid);
                s4 = (j < COLUMN_NUM-1) && (cells[i][j+1].type == liquid);
                s = s1 + s2 + s3 + s4;
                
                if(s1) cells[i][j].v -= d * (s1 / s);
                if(s2) cells[i+1][j].v += d * (s2 / s);
                if(s3) cells[i][j].u -= d * (s3 / s);
                if(s4) cells[i][j+1].u += d * (s4 / s);
                
                if(s) cells[i][j].p += (d / s) * (DENSITY * CELL_SIZE_M * TIME_STEP_HZ);
            }
}

void advection(SDL_Surface* surface, cell** cells){
    double v_bar;
    int x, y;
    cell** copy = (cell**)malloc(ROW_NUM * sizeof(cell*));
    
    for (int i = 0; i < ROW_NUM; ++i) {
        copy[i] = (cell*)malloc(COLUMN_NUM * sizeof(cell));
        for (int j = 0; j < COLUMN_NUM; ++j)
            copy[i][j] = cells[i][j];
    }
    
    for(int i = 0; i < ROW_NUM; ++i){
        for(int j = 0; j < COLUMN_NUM; ++j){
            if(copy[i][j].type == solid) continue;
            
//            if(i == 0 && j == COLUMN_NUM-1)
//                v_bar = cells[i][j].v;
//            else if(i == 0)
//                v_bar = (cells[i][j].v + cells[i][j+1].v) / 2;
//            else if(j == COLUMN_NUM-1)
//                v_bar = (cells[i][j].v + cells[i-1][j].v) / 2;
//            else
//                v_bar = (cells[i][j].v + cells[i][j+1].v + cells[i-1][j].v + cells[i-1][j+1].v) / 4;
//            y = (i * CELL_SIZE_M + v_bar / TIME_STEP_HZ) / CELL_SIZE_M;
//            x = (j * CELL_SIZE_M + copy[i][j].u / TIME_STEP_HZ) / CELL_SIZE_M;
//            y = (y < 0)? 0 : ((y >= ROW_NUM)? ROW_NUM - 1 : y);
//            x = (x < 0)? 0 : ((x >= COLUMN_NUM)? COLUMN_NUM - 1 : x);
//            
//            if(cells[i][j].type != solid)
//                cells[i][j] = copy[y][x];
            
            y = (i * CELL_SIZE_M + copy[i][j].v / TIME_STEP_HZ) / CELL_SIZE_M;
            x = (j * CELL_SIZE_M + copy[i][j].u / TIME_STEP_HZ) / CELL_SIZE_M;
            y = (y < 0)? 0 : ((y >= ROW_NUM)? ROW_NUM - 1 : y);
            x = (x < 0)? 0 : ((x >= COLUMN_NUM)? COLUMN_NUM - 1 : x);
            
            if(cells[y][x].type != solid) cells[y][x] = copy[i][j];
        }
    }
    free(copy);
}

void draw_grid(SDL_Surface* surface){
    for(int i = 0; i < COLUMN_NUM; ++i){
        SDL_Rect column = (SDL_Rect){i * CELL_SIZE, 0, LINE_WIDTH, WINDOW_HEIGHT};
        SDL_FillRect(surface, &column, GRID_COLOR);
    }
    for(int i = 0; i < ROW_NUM; ++i){
        SDL_Rect row = (SDL_Rect){0, i * CELL_SIZE, WINDOW_WIDTH, LINE_WIDTH};
        SDL_FillRect(surface, &row, GRID_COLOR);
    }
}

void draw_cell(SDL_Surface* surface, int x, int y, Uint32 color){
    SDL_Rect cell = (SDL_Rect){(x)*CELL_SIZE, (ROW_NUM-1-y)*CELL_SIZE, CELL_SIZE, CELL_SIZE};
    SDL_FillRect(surface, &cell, color);
}
void draw_cells(SDL_Surface* surface,  cell** cells){
    for(int i = 0; i < ROW_NUM; ++i){
        for(int j = 0; j < COLUMN_NUM; ++j){
            switch(cells[i][j].type){
                case liquid:
                    if(PRESSURE_MODE)
                        draw_cell(surface, j, i, cells[i][j].p);
                    else
                        draw_cell(surface, j, i, (cells[i][j].v * cells[i][j].v + cells[i][j].u * cells[i][j].u) );
                    break;
                case solid: draw_cell(surface, j, i, SOLID_COLOR); break;
                default: break;
            }
        }
    }
}

