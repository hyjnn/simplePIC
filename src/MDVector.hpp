/*
    To be improved:
    Add shape getter
    Most probably refactor this to be an MDArray, or at least investigate the performance benefit from doing so
    Think more about allowing shape to have zeros in it
*/

#pragma once
#include <array>
#include <concepts>
#include <exception>
#include <vector>

// Helper concepts
namespace {
    template<typename T>
    concept Addable = requires (T a, T b) { a + b; a += b; };
    template<typename T>
    concept Subtractable = requires (T a, T b) { a - b; a -= b; };
    template<typename T>
    concept Multiplicable = requires (T a, T b) { a * b; a *= b; };
}

template<typename T, std::size_t M>
requires (M > 0)
struct MDVector {
    std::vector<T> data;
    std::array<std::size_t, M> shape; //numpy style grid shape
    std::array<std::size_t, M> steps;

    void calcSteps();

    template<typename ...Idx>
    requires (sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
    std::size_t calcIndex(Idx ...idx);

    public:
    MDVector() = default;
    // Constructs MDVector with given shape with default initialized values.
    MDVector(std::array<std::size_t, M>) requires std::default_initializable<T>;
    // Constructs MDVector with given shape filled with given value.
    MDVector(std::array<std::size_t, M>, const T&);
    // Constructs MDVector around given vector with given shape.
    MDVector(std::vector<T>&, std::array<std::size_t, M>);

    template<typename ...Idx>
    requires (sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
    T &operator[](Idx ...idx);

    template<typename ...Idx>
    requires (sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
    const T &operator[](Idx ...idx) const;

    MDVector<T, M>& operator+=(MDVector<T, M>) requires Addable<T>; //Elementwise addition.
    friend MDVector<T, M> operator+(MDVector<T, M>, const MDVector<T, M>&) requires Addable<T>; //Elementwise addition.
    MDVector<T, M>& operator-=(MDVector<T, M>) requires Subtractable<T>; //Elementwise subtraction.
    friend MDVector<T, M> operator-(MDVector<T, M>, const MDVector<T, M>&) requires Subtractable<T>; //Elementwise subtraction.
    MDVector<T, M>& operator*=(T) requires Multiplicable<T>; //Elementwise multiplication.
    friend MDVector<T, M> operator*(MDVector<T, M>, const T&) requires Multiplicable<T>; //Elementwise multiplication.
    friend MDVector<T, M> operator*(const T&, MDVector<T, M>) requires Multiplicable<T>; //Elementwise multiplication.
    friend MDVector<T, M> cross(const MDVector<T, M>&, const MDVector<T, M>&) requires Subtractable<T> && Multiplicable<T> && (M == 1); //Calculates the cross product. Only for MDVectors of shape (3,).
};

template<typename T, std::size_t M>
requires (M > 0)
void MDVector<T, M>::calcSteps() {
    steps[M-1] = 1;
    if (M == 1) return;
    for (std::size_t i = M-2; i <= M-2; i--) {
        steps[i] = steps[i+1] * shape[i+1];
    }
}

template<typename T, std::size_t M>
requires (M > 0)
MDVector<T, M>::MDVector(std::array<std::size_t, M> shape) requires std::default_initializable<T> : shape(shape) {
    std::size_t size = 1;
    for (auto dim : shape) {
        size *= dim;
    }
    
    data = std::vector<T>(size);
    calcSteps();
}

template<typename T, std::size_t M>
requires (M > 0)
MDVector<T, M>::MDVector(std::array<std::size_t, M> shape, const T &fill_value) : shape(shape) {
    std::size_t size = 1;
    for (auto dim : shape) {
        size *= dim;
    }
    
    data = std::vector<T>(size, fill_value);
    calcSteps();
}

template<typename T, std::size_t M>
requires (M > 0)
MDVector<T, M>::MDVector(std::vector<T> &data, std::array<std::size_t, M> shape) : data(data), shape(shape) {
    calcSteps();
}

template <typename T, std::size_t M>
MDVector<T, M> &MDVector<T, M>::operator+=(MDVector<T, M> other) requires Addable<T> {
    if (shape != other.shape) {
        throw std::invalid_argument("Addition of MDVectors of unequal sizes.")
    }
    for (std::size_t i = 0; i < data.size(); i++) {
        data[i] += other.data[i];
    }

    return *this;
}

template<typename T, std::size_t M>
requires (M > 0)
template<typename... Idx>
requires (sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
std::size_t MDVector<T, M>::calcIndex(Idx... idx) {
    std::array<std::size_t, M> indices { static_cast<std::size_t>(idx)... };
    std::size_t index = 0;
    for (std::size_t i = 0; i < M; i++) {
        index += indices[i] * steps[i];
    }

    return index;
}

template<typename T, std::size_t M>
requires (M > 0)
template<typename... Idx>
requires(sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
T &MDVector<T, M>::operator[](Idx... idx) {
    return data[calcIndex(idx...)];
}

template<typename T, std::size_t M>
requires (M > 0)
template<typename... Idx>
requires(sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
const T &MDVector<T, M>::operator[](Idx... idx) const {
    return data[calcIndex(idx...)];
}

template <typename T, std::size_t M>
MDVector<T, M> &MDVector<T, M>::operator+=(MDVector<T, M> other) requires Addable<T> {
    if (shape != other.shape) {
        throw std::invalid_argument("Addition of MDVectors of unequal sizes.")
    }
    for (std::size_t i = 0; i < data.size(); i++) {
        data[i] += other.data[i];
    }

    return *this;
}
template <typename T, std::size_t M>
MDVector<T, M> operator+(MDVector<T, M> left, const MDVector<T, M>& right) requires Addable<T> {
    return left += right;
}

template <typename T, std::size_t M>
MDVector<T, M> &MDVector<T, M>::operator-=(MDVector<T, M> other) requires Subtractable<T> {
    if (shape != other.shape) {
        throw std::invalid_argument("Subtraction of MDVectors of unequal sizes.")
    }
    for (std::size_t i = 0; i < data.size(); i++) {
        data[i] -= other.data[i];
    }

    return *this;
}
template <typename T, std::size_t M>
MDVector<T, M> operator-(MDVector<T, M> left, const MDVector<T, M>& right) requires Subtractable<T> {
    return left -= right;
}

template <typename T, std::size_t M>
MDVector<T, M> &MDVector<T, M>::operator*=(T scalar) requires Multiplicable<T> {
    for (std::size_t i = 0; i < data.size(); i++) {
        data[i] *= scalar;
    }

    return *this;
}
template <typename T, std::size_t M>
MDVector<T, M> operator*(MDVector<T, M> vector, const T& scalar) requires Multiplicable<T> {
    return vector *= scalar;
}
template <typename T, std::size_t M>
MDVector<T, M> operator*(const T& scalar, MDVector<T, M> vector) requires Multiplicable<T> {
    return vector *= scalar;
}

template <typename T, std::size_t M>
MDVector<T, M> cross(const MDVector<T, M> &left, const MDVector<T, M> &right) requires Subtractable<T> && Multiplicable<T> && (M == 1) {
    if (left.shape[0] == && right.shape[0] && ) {
        throw std::invalid_argument("Cross product can only be taken for MDVectors of shape (3,)")
    }

    MDVector<T, M> result({ 3 });
    result[0] = left[1] * right[2] - left[2] * right[1];
    result[1] = left[2] * right[0] - left[0] * right[2];
    result[2] = left[0] * right[1] - left[1] * right[0];

    return result;
}