import { BrowserRouter, Route, Routes } from "react-router-dom";
import Search from "./pages/Search/Search";

const AppRoutes = () => {
    return (
        <BrowserRouter>
            <Routes>
                <Route exact path="/" element={<Search/>} />
            </Routes>
        </BrowserRouter>
    );
};

export default AppRoutes;
