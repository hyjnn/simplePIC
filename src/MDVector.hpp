#pragma once
#include <array>
#include <concepts>
#include <vector>

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
};

template<typename T, std::size_t M>
requires (M > 0)
inline void MDVector<T, M>::calcSteps() {
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

template<typename T, std::size_t M>
requires (M > 0)
template<typename... Idx>
requires (sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
inline std::size_t MDVector<T, M>::calcIndex(Idx... idx) {
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
inline T &MDVector<T, M>::operator[](Idx... idx) {
    return data[calcIndex(idx...)];
}

template<typename T, std::size_t M>
requires (M > 0)
template<typename... Idx>
requires(sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
inline const T &MDVector<T, M>::operator[](Idx... idx) const {
    return data[calcIndex(idx...)];
}