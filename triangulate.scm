(import srfi-133)

;; vectors
(import srfi-1)

;; filter
(import (chicken io)
        (chicken random)
        (chicken port)
        (chicken process-context)
        (chicken string)
        (chicken sort)
        (chicken format))
(import (prefix sdl2 "sdl2:"))

(define pixels-per-point 1)
(define epsilon 0.001)
(define *color-back* (sdl2:make-color 8 26 0))
(define *color-point* (sdl2:make-color 230 238 255))
(define *color-circle* (sdl2:make-color 179 189 0))
(define *color-rect* (sdl2:make-color 255 77 77))
(define *color-triangle* (sdl2:make-color 10 230 20))
(define-record-type trig
 (make-trig p1 p2 p3 done?)
 trig?
 (p1 trig-p1 trig-p1-set!)
 (p2 trig-p2 trig-p2-set!)
 (p3 trig-p3 trig-p3-set!)
 (done? trig-done? trig-done-set!))

(define-record-type circle
 (make-circle pc rad)
 circle?
 (pc circle-pc circle-pc-set!)
 (rad circle-rad circle-rad-set!))

(define-record-printer (trig t out)
 (fprintf out
          "#trig(done?:~s): ((~s,~s),(~s,~s),(~s,~s))"
          (trig-done? t)
          (vector-ref (trig-p1 t) 0)
          (vector-ref (trig-p1 t) 1)
          (vector-ref (trig-p2 t) 0)
          (vector-ref (trig-p2 t) 1)
          (vector-ref (trig-p3 t) 0)
          (vector-ref (trig-p3 t) 1)))

(define-record-printer (circle c out)
 (fprintf out
          "#circ: ((~s,~s)/~s))"
          (vector-ref (circle-pc c) 0)
          (vector-ref (circle-pc c) 1)
          (circle-rad c)))

(define (has-shared-verts? a b)
 (let ((a1 (trig-p1 a))
       (a2 (trig-p2 a))
       (a3 (trig-p3 a))
       (b1 (trig-p1 b))
       (b2 (trig-p2 b))
       (b3 (trig-p3 b)))
  (or (v= a1 b1)
      (v= a1 b2)
      (v= a1 b3)
      (v= a2 b1)
      (v= a2 b2)
      (v= a2 b3)
      (v= a3 b1)
      (v= a3 b2)
      (v= a3 b3))))

(define-record-type edge
 (make-edge p1 p2)
 edge?
 (p1 edge-p1 edge-p1-set!)
 (p2 edge-p2 edge-p2-set!))

(define-record-printer (edge e out)
 (fprintf out "#edge: ~s <-> ~s" (edge-p1 e) (edge-p2 e)))

;; HELPERS
(define (edge= e1 e2)
 (or (and (v= (edge-p1 e1) (edge-p1 e2))
          (v= (edge-p2 e1) (edge-p2 e2)))
     (and (v= (edge-p1 e1) (edge-p2 e2))
          (v= (edge-p2 e1) (edge-p1 e2)))))

(define (square x)
 (* x x))

(define (triangle-to-edges t)
 (list (make-edge (trig-p1 t) (trig-p2 t))
       (make-edge (trig-p2 t) (trig-p3 t))
       (make-edge (trig-p3 t) (trig-p1 t))))

(define (v= a b)
 (vector= = a b))

;; OBJ FILE ROUTINES 
;; reads one line of an obj file and adds the point to vector of points
(define (process-obj-line l verts)
 (cond
  ((and (substring=? "v " l)
        (> (string-length l) 3))
   (let ((fields (string-split l " ")))
    (if (> (length fields) 1)
     (vector-append
      verts
      (vector (list->vector (map (lambda (x)
                                  (+ 10 (* 10 (string->number x))))
                                 (cdr fields)))))
     verts)))
  (else
   verts)))

;; read an obj file and collect the vertexes
(define (read-obj-file file)
 (let ((fh (open-input-file file))
       (verts '#()))
  (let loop ((c (read-line fh)))
   (if (eof-object? c)
    (begin
     (close-input-port fh)
     verts)
    (begin
     (set! verts (process-obj-line c verts))
     (loop (read-line fh)))))
  (format #f "verts :~A~%" verts)
  verts))

(define reader-test-files (vector "data/airboat.obj" "data/cube.obj"))
(define (reader-tests)
 (vector-map read-obj-file (vector (vector-ref reader-test-files 1))))

;; TRIANGULATOR ROUTINES
;; returns a triangle  that totally surround the point cloud
;; the fudge factor should be a nonnegative integer
(define (get-bounding-triangle pc fudge)
 (let* ((mmax (get-min-maxs pc))
        (d1 (car mmax))
        (d2 (cadr mmax))
        (xmid (ceiling (inexact->exact (/ (+ (car d1) (cdr d1)) 2))))
        (dx (- (car d1) (cdr d1)))
        (ymid (ceiling (inexact->exact (/ (+ (car d2) (cdr d2)) 2))))
        (dy (- (car d2) (cdr d2)))
        (dmax (max dx dy)))
  (make-trig (vector (ceiling (- xmid (* fudge dmax))) (- ymid dmax))
             (vector xmid (ceiling (+ ymid (* fudge dmax))))
             (vector (ceiling (+ xmid (* fudge dmax))) (- ymid dmax))
             #f)))

(define (triangulate pc)
 (let* ((supertriangle (get-bounding-triangle pc 1.35))
        ; start the trig list with the supertriangle")
        (triglist (list supertriangle))
        ; for each point p
        (pointiter (lambda (acc p)
		    ; init empty edge buffer
                    (let* ((edgebuf '())
                           ; for every triangle (accumulate also some couners)
                           (trigsprocessed (fold
                                            (lambda (t acc1)
					     ; find the circumcircle
                                             (let ((ccircle (circumcircle t))
                                                   (slen (length triglist)))
                                              ; save the length of the triglist to find what we removed
                                              (if (in-circle? p ccircle)
                                               ; if the point p is in the circumcircle
                                               (begin
                                                (set! edgebuf
                                                 (append edgebuf
                                                         (triangle-to-edges t)))
                                                ; add the triangle's edges to edgebuf
                                                (set! triglist
                                                 (filter (lambda (x)
                                                          (not (eq? x t)))
                                                         ; remove this triangle
                                                         triglist))))
                                              ; return a pair of the total trigs processed and how many were removed
                                              (cons (+ (car acc1) 1)
                                                    (- slen (length triglist)))))
                                            (cons 0 0)
                                            ; just the initial accumulator for the analytics
                                            triglist)))
                     ;; delete all doubly specified edges from edgebuf, which will only leave the enclosing polygon
                     (set! edgebuf (dedupe edgebuf edge=))
                     ;; add to the triangle list the triangles formed by the point and the edges of the polygon
                     (set! triglist
                      (append triglist
                              (map (lambda (x)
                                    (make-trig (edge-p1 x) (edge-p2 x) p #f))
                                   edgebuf)))))))
  (vector-fold pointiter 0 pc)
  ;; remove any triangles that share any supertriangle vertices
  (filter (lambda (x)
           (not (has-shared-verts? x supertriangle)))
          triglist)))
  ;; no need to remove the supertriangle vertices from the vertex list since we never added

;; pass a list and an equality function
(define (dedupe e eqf)
 (if (null? e)
  '()
  (cons (car e)
        (dedupe (filter (lambda (x)
                         (not (eqf x (car e))))
                        (cdr e))
                eqf))))

;; circumcircle takes a point and a triangle and returns
;; a circle that passes from the three points of the rectangle
(define (circumcircle trig)
 (format #t "doing trig :~A~%" trig)
 (let* ((p1 (trig-p1 trig))
        (p2 (trig-p2 trig))
        (p3 (trig-p3 trig))
        (x1 (vector-ref p1 0))
        (y1 (vector-ref p1 1))
        (x2 (vector-ref p2 0))
        (y2 (vector-ref p2 1))
        (x3 (vector-ref p3 0))
        (y3 (vector-ref p3 1))
        (dyp1p2 (abs (- y1 y2)))
        (dyp2p3 (abs (- y2 y3)))
        (xc 0)
        (yc 0)
        (r 0))
  (cond
   ((and (< dyp1p2 epsilon)
         (< dyp2p3 epsilon))
    'COLINEAR_POINTS)
   ((< dyp1p2 epsilon)
    (let ((m2 (/ (- (- x3 x2)) (- y3 y2)))
          (mx2 (/ (+ x2 x3) 2))
          (my2 (/ (+ y2 y3) 2)))
     (set! xc (/ (+ x2 x1) 2))
     (set! yc (+ (* m2 (- xc mx2)) my2))))
   ((< dyp2p3 epsilon)
    (let ((m1 (/ (- (- x2 x1)) (- y2 y1)))
          (mx1 (/ (+ x1 x2) 2))
          (my1 (/ (+ y1 y2) 2)))
     (set! xc (/ (+ x3 x2) 2))
     (set! yc (+ (* m1 (- xc mx1)) my1))))
   (else
    (let ((m1 (/ (- (- x2 x1)) (- y2 y1)))
          (mx1 (/ (+ x1 x2) 2))
          (my1 (/ (+ y1 y2) 2))
          (m2 (/ (- (- x3 x2)) (- y3 y2)))
          (mx2 (/ (+ x2 x3) 2))
          (my2 (/ (+ y2 y3) 2)))
     (if (> (- m1 m2) epsilon)
      ; guard from div zero
      (set! xc (/ (+ (- (* m1 mx1) (* m2 mx2)) my2 (- my1)) (- m1 m2)))
      (set! xc (/ (+ (- (* m1 mx1) (* m2 mx2)) my2 (- my1)) epsilon)))
     (if (> dyp1p2 dyp2p3)
      (set! yc (+ (* m1 (- xc mx1)) my1))
      (set! yc (+ (* m2 (- xc mx2)) my2))))))
  (make-circle (vector xc yc) (sqrt (+ (square (- x2 xc)) (square (- y2 yc)))))))

(define (in-circle? p c)
 (let* ((xp (vector-ref p 0))
        (yp (vector-ref p 1))
        (pc (circle-pc c))
        (r (circle-rad c))
        (xc (vector-ref pc 0))
        (yc (vector-ref pc 1)))
  (<= (- (sqrt (+ (square (- xp xc)) (square (- yp yc)))) r) epsilon)))

;; POINTCLOUD ROUTINES
;; sort by axis
(define (find-sizes pc)
 (let* ((mmax (get-min-maxs pc))
        (d1 (car mmax))
        (d2 (cadr mmax))
        (d3 (caddr mmax)))
  (vector (- (car d1) (cdr d1)) (- (car d2) (cdr d2)) (- (car d3) (cdr d3)))))

(define (sort-by-x vec)
 (sort vec
       (lambda (a b)
        (< (vector-ref a 0) (vector-ref b 0)))))

(define (sort-by-y vec)
 (sort vec
       (lambda (a b)
        (< (vector-ref a 1) (vector-ref b 1)))))

(define (sort-by-z vec)
 (sort vec
       (lambda (a b)
        (< (vector-ref a 2) (vector-ref b 2)))))

;; returns a list with pairs of min max in each dim
(define (get-min-maxs pc)
 (let ((pcx (sort-by-x pc))
       (pcy (sort-by-y pc))
       (pcz (sort-by-z pc))
       (last (- (vector-length pc) 1)))
  (list
   (cons (vector-ref (vector-ref pcx last) 0)
         (vector-ref (vector-ref pcx 0) 0))
   (cons (vector-ref (vector-ref pcy last) 1)
         (vector-ref (vector-ref pcy 0) 1))
   (cons (vector-ref (vector-ref pcz last) 2)
         (vector-ref (vector-ref pcz 0) 2)))))

;; return a pair of the minimum and the maximum number in all dimensions
(define (get-global-min-max pc)
 (let* ((mmax (get-min-maxs pc))
        (dmax (vector-fold (lambda (acc v)
                            (max v acc))
                           0
                           mmax))
        (dmin (vector-fold (lambda (acc v)
                            (min v acc))
                           0
                           mmax)))
  ; min dimension size
  (cons dmax dmin)))

;; random poingcloud generator for testing
(define (gen-rand-points n maxsz)
 (vector-unfold (lambda (i)
                 (vector (+ (/ maxsz 2) (pseudo-random-integer maxsz))
                         (+ (/ maxsz 2) (pseudo-random-integer maxsz))
                         (+ (/ maxsz 2) (pseudo-random-integer maxsz))))
                n))

;; GRAPHICS methods
(define (draw-points vec)
 (let* ((sizes (find-sizes vec))
        (screensz (vector-map (lambda (a)
                               (* 3 (* pixels-per-point a)))
                              sizes))
        (window (sdl2:create-window! "chicken!triangulate!"
                                     0
                                     0
                                     (ceiling (vector-ref screensz 0))
                                     (ceiling (vector-ref screensz 1))))
        (renderer (sdl2:create-renderer! window -1 '(accelerated)))
        (points (vector-map
                 (lambda (p)
                  (sdl2:make-point (ceiling (vector-ref p 0))
                                   (ceiling (vector-ref p 1))))
                 vec))
        (alltrigs (triangulate vec)))
  (sdl2:set-main-ready!)
  (sdl2:init! '(video))
  ;; background
  (sdl2:render-draw-color-set! renderer *color-back*)
  (sdl2:render-draw-blend-mode-set! renderer 'none)
  (sdl2:fill-rect! (sdl2:window-surface window) #f *color-back*)
  ;; points
  (sdl2:render-draw-color-set! renderer *color-point*)
  (sdl2:render-draw-points! renderer points)
  ;; little triangles
  (map (lambda (x)
        (draw-triangle renderer x))
       alltrigs)))
  ;(sdl2:delay! 13000)
  ;(sdl2:quit!)))

(define (draw-triangle renderer t)
 (sdl2:render-draw-color-set! renderer *color-triangle*)
 (sdl2:render-draw-lines!
  renderer
  (map (lambda (p)
        (sdl2:make-point (vector-ref p 0) (vector-ref p 1)))
       (list (trig-p1 t) (trig-p2 t) (trig-p3 t))))
 (sdl2:render-present! renderer))
