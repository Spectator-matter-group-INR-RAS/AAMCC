//
// Created by Artem Novikov on 17.02.2024.
//

#ifndef HANDLER_MYFERMIBREAKUP_UTILITIES_SINGLETON_H_
#define HANDLER_MYFERMIBREAKUP_UTILITIES_SINGLETON_H_

#include <memory>

template <typename T>
class Singleton {
 public:
  Singleton() {
#if __cplusplus >= 202002L
    if (instance_ == nullptr) [[unlikely]] {
#else
    if (instance_ == nullptr) {
#endif
      instance_ = std::make_unique<T>();
    }
  }

  template <typename... Args>
  Singleton(Args&& ... args) {
    instance_ = std::move(std::make_unique<T>(std::forward<Args>(args)...));
  }

  Singleton(T* ptr) {
    instance_.reset(ptr);
  }

  ~Singleton() = default;

  template <typename... Args>
  static void Reset(Args&& ... args) {
    instance_ = std::make_unique<T>(std::forward<Args>(args)...);
  }

  static void Reset(T* ptr) {
    instance_.reset(ptr);
  }

  static const T& Instance() {
#if __cplusplus >= 202002L
    if (instance_ == nullptr) [[unlikely]] {
#else
    if (instance_ == nullptr) {
#endif
      instance_ = std::make_unique<T>();
    }
    return *instance_;
  }

  T& operator*() {
    return *instance_;
  }

  const T& operator*() const {
    return *instance_;
  }

  T* operator->() {
    return instance_.get();
  }

  const T& operator->() const {
    return *instance_;
  }

 private:
  static inline std::unique_ptr<T> instance_;
};

#endif //HANDLER_MYFERMIBREAKUP_UTILITIES_SINGLETON_H_
