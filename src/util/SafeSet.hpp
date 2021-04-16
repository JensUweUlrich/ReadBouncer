#pragma once

#include <condition_variable>
#include <mutex>
#include <set>
#include <iterator>
#include <type_traits>

#ifndef SAFESET_HPP_
#define SAFESET_HPP_

template < class T,
           class Enabled = std::enable_if_t< std::is_move_assignable< T >::value && //
                                             std::is_move_constructible< T >::value > >
class SafeSet
{
    private:
        std::set< T > s;
        std::mutex      m;
        // Very high max_size - default to std::numeric_limits< size_t >::max() - virtually no limit on size of the map
        size_t                  max_size;
        bool                    push_over = false;
        std::condition_variable cv_push;
        std::condition_variable cv_pop;

    public:
        SafeSet()
        : max_size( std::numeric_limits< size_t >::max() )
        {
        }


        SafeSet( size_t max )
        : max_size( max )
        {
        }

        void set_max_size( size_t max )
        {
            std::lock_guard< std::mutex > lock( m );
            max_size = max;
        }

    

        void insert( T t )
        {
            // accquire mutex to modify queue
            std::unique_lock< std::mutex > lock( m );
            // insert pair in the map
            s.insert(t);
            // notify other thread that something is in the map
            cv_pop.notify_one();
        }

        void erase( T t)
        {
            std::unique_lock< std::mutex > lock( m );
            s.erase(t);
            cv_push.notify_one();
        }

        void erase(typename std::set<T>::iterator t)
        {
            std::unique_lock< std::mutex > lock(m);
            s.erase(t);
            cv_push.notify_one();
        }

        typename std::set<T>::iterator end()
        {
            std::lock_guard< std::mutex > lock(m);
            return s.end();
        }

        void notify_push_over()
        {
            std::lock_guard< std::mutex > lock( m );
            push_over = true;
            cv_pop.notify_all();
        }

        int size()
        {
            std::lock_guard< std::mutex > lock( m );
            return s.size();
        }

        bool empty()
        {
            std::lock_guard< std::mutex > lock( m );
            return s.empty();
        }

        typename std::set<T>::iterator find(T t)
        {
            std::lock_guard< std::mutex > lock(m);
            return s.find(t);
        }

        bool contains(T t)
        {
            std::lock_guard< std::mutex > lock(m);
            std::set<T>::iterator it = s.find(t);
            if (it != s.end())
                return true;
            else
                return false;
        }
};

#endif // SAFESET_HPP_