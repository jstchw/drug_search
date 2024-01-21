import React, { useContext } from "react";
import { SearchHistoryContext } from "../../contexts/SearchHistoryContext";
import {Button, ListGroup, OverlayTrigger, Popover} from "react-bootstrap";
import { Trash, Binoculars, GenderIntersex, Globe, Baby } from "@phosphor-icons/react"
import { mapParamToLabel } from "../../utils/utils";
import {searchSex, searchTypes} from "../../constants";
import {useSearch} from "../../hooks/useSearch";
import {URLParams} from "../../types";

const SearchHistory = ({ setShowOptionModal }: { setShowOptionModal: React.Dispatch<boolean>}) => {
    const context = useContext(SearchHistoryContext)
    const { executeSearch } = useSearch('URLParams')

    if (!context) {
        return null
    } else if (context.searchHistory.length === 0) {
        return (
            <>
                <div>
                    <ListGroup>
                        <ListGroup.Item>
                            <div className={'fs-5'} style={{width: 'fit-content'}}>No search history</div>
                        </ListGroup.Item>
                    </ListGroup>
                </div>
            </>
        )
    }

    const { searchHistory, clearSearchHistory } = context;

    const handleRedirect = (terms: string[], searchOptions: URLParams) => {
        executeSearch(terms,
            {
                terms: [],
                searchBy: searchOptions.searchBy,
                sex: searchOptions.sex,
                age: searchOptions.age,
                country: searchOptions.country
            })
        setShowOptionModal(false)
    }

    return (
        <>
            <div>
                <ListGroup>
                {searchHistory.map((item, index) => {
                    const popover = (
                        <Popover>
                            <Popover.Header style={{backgroundColor: 'transparent'}}>
                                <div className={'d-flex align-items-center fs-5'}>
                                    <Binoculars weight={'light'}/>
                                    <div className={'vr mx-2'}/>
                                    <span>{item.searchBy && mapParamToLabel(item.searchBy, searchTypes)}</span>
                                </div>
                            </Popover.Header>
                            <Popover.Body>
                                <div className={'d-flex align-items-center fs-5'}>
                                    <GenderIntersex weight={'light'}/>
                                    <div className={'vr mx-2'}/>
                                    <span>{item.sex ? mapParamToLabel(item.sex, searchSex) : <span className={'text-secondary'}>Not specified</span>}</span>
                                </div>
                                <div className={'d-flex align-items-center fs-5'}>
                                    <Baby weight={'light'}/>
                                    <div className={'vr mx-2'}/>
                                    <span>{item.age && item.age.min ? item.age.min :
                                        <span className={'text-secondary'}>N/a</span>}</span>
                                    <span> - </span>
                                    <span>{item.age && item.age.max ? item.age.max : <span className={'text-secondary'}>N/a</span>}</span>
                                </div>
                                <div className={'d-flex align-items-center fs-5'}>
                                    <Globe weight={'light'}/>
                                    <div className={'vr mx-2'}/>
                                    <span>{item.country ? item.country : <span className={'text-secondary'}>Not specified</span>}</span>
                                </div>
                            </Popover.Body>
                        </Popover>
                    );

                    return (
                        <ListGroup.Item key={index}>
                            <OverlayTrigger overlay={popover} key={index} trigger={['hover', 'focus']} placement="bottom">
                                <div onClick={() => handleRedirect(item.terms, item)} className={'fs-5'} key={index} style={{width: 'fit-content', cursor: 'pointer'}}>
                                    {item.terms.map(term => term.charAt(0).toUpperCase() + term.slice(1)).join(' & ')}
                                </div>
                            </OverlayTrigger>
                        </ListGroup.Item>
                    );
                })}
                    <ListGroup.Item className={'d-flex justify-content-end'}>
                        <Button variant={'outline-primary'} className={'d-flex align-items-center'} onClick={clearSearchHistory}>
                            <Trash/>
                        </Button>
                    </ListGroup.Item>
                </ListGroup>
            </div>
        </>
    )
}

export default SearchHistory;