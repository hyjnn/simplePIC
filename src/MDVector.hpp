#pragma once
#include <array>
#include <concepts>
#include <vector>

template<typename T, std::size_t M>
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
    MDVector(std::array<std::size_t, M>) requires std::default_initializable<T>;
    MDVector(std::vector<T>&, std::array<std::size_t, M>);

    template<typename ...Idx>
    requires (sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
    T &operator[](Idx ...idx);

    template<typename ...Idx>
    requires (sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
    const T &operator[](Idx ...idx) const;
};

template<typename T, std::size_t M>
inline void MDVector<T, M>::calcSteps() {
    steps[M-1] = 1;
    for (std::size_t i = M-2; i >= 0; i--) {
        steps[i] = steps[i+1] * shape[i+1];
    }
}

template<typename T, std::size_t M>
MDVector<T, M>::MDVector(std::array<std::size_t, M> shape) requires std::default_initializable<T> : shape(shape) {
    std::size_t size = 1;
    for (auto dim : shape) {
        size *= dim;
    }
    
    data = std::vector<T>(size);
}

template<typename T, std::size_t M>
MDVector<T, M>::MDVector(std::vector<T> &data, std::array<std::size_t, M> shape) : data(data), shape(shape) {
    calcSteps();
}

template<typename T, std::size_t M>
template<typename... Idx>
requires (sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
inline std::size_t MDVector<T, M>::calcIndex(Idx... idx) {
    indices = std::array<std::size_t, M>{ std::static_cast<std::size_t>(idx...) }
    index = 0;
    for (std::size_t i = 0; i < M; i++) {
        index += indices[i] * steps[i];
    }

    return index;
}

template<typename T, std::size_t M>
template<typename... Idx>
requires(sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
inline T &MDVector<T, M>::operator[](Idx... idx) {
    return data[calcIndex(idx...)];
}

template<typename T, std::size_t M>
template<typename... Idx>
requires(sizeof...(Idx) == M) && (std::is_convertible_v<Idx, std::size_t> && ...)
inline const T &MDVector<T, M>::operator[](Idx... idx) const {
    return data[calcIndex(idx...)];
}